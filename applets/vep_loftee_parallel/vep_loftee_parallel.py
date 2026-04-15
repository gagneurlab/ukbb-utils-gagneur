#!/usr/bin/env python3

"""
TO BUILD: dx build applets/vep_loftee_parallel/ --destination <PROJECT-ID>:<PATH-TO-APPLET-DEST>
TO BUILD: dx build applets/vep_loftee_parallel/ --destination project-REDACTED:/users/sl/ -f
"""

import os
import subprocess
import dxpy
import logging
import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Single source of truth for paths
WORK_DIR = "/home/dnanexus"
VEP_DATA = os.path.join(WORK_DIR, "vep_data")

def run(cmd):
    logger.info(f"$ {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def parquet_to_vcf(parquet_path, vcf_path):
    """Converts input parquet to a minimal VCF for VEP."""
    logger.info(f"Converting {parquet_path} to VCF...")
    df = pl.read_parquet(parquet_path)

    # 1. Apply your local Karyotypic Sorting logic!
    df = df.with_columns(
        chrom_sort_key=pl.when(
            pl.col("chrom").str.replace("(?i)^chr", "").str.to_uppercase() == "X"
        ).then(23)
        .when(
            pl.col("chrom").str.replace("(?i)^chr", "").str.to_uppercase() == "Y"
        ).then(24)
        .when(
            pl.col("chrom").str.replace("(?i)^chr", "").str.to_uppercase() == "MT"
        ).then(25)
        .otherwise(
            # strict=False safely handles unmapped contigs (like GL000191) by making them Null
            pl.col("chrom").str.replace("(?i)^chr", "").cast(pl.Int64, strict=False)
        )
    ).sort(["chrom_sort_key", "pos"])

    # 2. Format columns
    df_vcf = df.select([
        pl.col("chrom").cast(pl.Utf8).alias("#CHROM"),
        pl.col("pos").alias("POS"),
        pl.col("id").fill_null(".").alias("ID"),
        pl.col("ref").alias("REF"),
        pl.col("alt").alias("ALT"),
        pl.lit(".").alias("QUAL"),
        pl.lit(".").alias("FILTER"),
        pl.lit(".").alias("INFO")
    ])

    # 3. Fast native write
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(df_vcf.write_csv(separator="\t", include_header=True))

def _parse_vep_extra_col(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Extracts strict LOFTEE fields from the VEP 'extra' column natively in Polars.
    This replaces the slow Python string-split loop to ensure the pipeline remains lazy.
    """
    logger.info("Extracting LOFTEE and standard fields from 'extra' column...")
    
    # Check if 'extra' exists in the lazy schema
    if "extra" not in lf.collect_schema().names():
        return lf

    extra_exprs = [
        # The (?:^|;) regex ensures we match the key exactly, whether it's at the start or middle
        pl.col("extra").str.extract(r"(?:^|;)LoF=([^;]*)", 1).alias("lof"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_filter=([^;]*)", 1).alias("lof_filter"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_flags=([^;]*)", 1).alias("lof_flags"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_info=([^;]*)", 1).alias("lof_info"),
        pl.col("extra").str.extract(r"(?:^|;)IMPACT=([^;]*)", 1).alias("impact")
    ]
    
    # Apply the extractions and drop the heavy original text column
    return lf.with_columns(extra_exprs).drop("extra")

def post_process_vep(vep_lf, metadata_parquet_path, gene_list=None, biotypes_filter=None, use_loftee=True):
    """
    Takes raw VEP TSV data, joins original variant metadata, applies filters,
    parses extras, generates dummies, and returns a processed LazyFrame.
    """
    logger.info("Loading original variant metadata for joining...")
    meta_lf = pl.scan_parquet(metadata_parquet_path)
    
    # Standardize VEP's variation column to match the Parquet's "id" column
    vep_lf = vep_lf.rename({"#Uploaded_variation": "id"})

    # 1. The Master Join (Restoring CHROM, POS, REF, ALT)
    logger.info("Joining VEP annotations with original variant metadata...")
    lf = meta_lf.join(vep_lf, on="id", how="inner")
    
    # Standardize all column names to lowercase
    lf = lf.rename({col: col.lower() for col in lf.collect_schema().names()})

    # 2. Extract specific fields from the 'extra' column
    lf = _parse_vep_extra_col(lf)

    # 3. Apply Gene & Biotype Filters
    if gene_list:
        logger.info(f"Filtering for {len(gene_list)} target genes...")
        lf = lf.filter(pl.col("gene").is_in(gene_list))

    if biotypes_filter:
        logger.info(f"Filtering biotypes: {biotypes_filter}")
        lf = lf.filter(pl.col("biotype").is_in(biotypes_filter))

    # 4. LOFTEE Processing 
    if use_loftee and "lof" in lf.collect_schema().names():
        logger.info("Processing LOFTEE dummy annotations...")
        lf = lf.with_columns([
            pl.col("lof").eq("HC").cast(pl.Int8).fill_null(0).alias("loftee_hc"),
            pl.col("lof").eq("LC").cast(pl.Int8).fill_null(0).alias("loftee_lc"),
            pl.col("lof").is_null().cast(pl.Int8).alias("loftee_hc_is_na"),
            pl.col("lof").is_null().cast(pl.Int8).alias("loftee_lc_is_na")
        ])

    # 5. Handle 'consequence' Dummies
    logger.info("One-hot encoding consequences...")
    lf = lf.with_row_index("row_nr")

    dummies = (
        lf.select(["row_nr", "consequence"])
        .with_columns(pl.col("consequence").str.split(","))
        .explode("consequence")
        .drop_nulls("consequence")
        .filter(pl.col("consequence") != "")
        .collect()
        .to_dummies(columns="consequence")
        .group_by("row_nr").max()
    )

    dummies = dummies.rename({c: c.lower() for c in dummies.columns})
    dummy_cols = [c for c in dummies.columns if c != "row_nr"]
    dummies = dummies.with_columns([pl.col(c).cast(pl.Int8) for c in dummy_cols])

    # 6. Final Deduplication
    processed_lf = (
        lf.join(dummies.lazy(), on="row_nr", how="left")
        .drop("row_nr")
        .unique()
        .rename({'gene': 'region'})
    )
    return processed_lf


@dxpy.entry_point("process_chunk")
def process_chunk(chunk_file, vep_version, use_loftee):
    # 1. DEFINE & INITIALIZE PATHS
    PLUGINS_DIR = os.path.join(VEP_DATA, "Plugins")
    
    logger.info(f"Initializing directories at {VEP_DATA}...")
    os.makedirs(VEP_DATA, exist_ok=True)
    os.makedirs(PLUGINS_DIR, exist_ok=True)

    # 2. PREPARE INPUT
    local_parquet = os.path.join(WORK_DIR, "input.parquet")
    logger.info(f"Downloading chunk...")
    dxpy.download_dxfile(chunk_file, local_parquet)

    # Place input VCF inside VEP_DATA so it's accessible via the cache mount
    local_vcf = os.path.join(VEP_DATA, "input.vcf")
    parquet_to_vcf(local_parquet, local_vcf)

    # DEBUG: Confirm the VCF is exactly where Docker expects it
    run(f"ls -lh {local_vcf}")

    # 3. DOWNLOAD & UNPACK CACHE
    if not os.path.exists(os.path.join(VEP_DATA, "homo_sapiens")):
        logger.info(f"Downloading VEP {vep_version} Cache from Ensembl FTP...")
        ftp_url = f"http://ftp.ensembl.org/pub/release-{vep_version}/variation/vep/homo_sapiens_vep_{vep_version}_GRCh38.tar.gz"
        cache_tar = os.path.join(WORK_DIR, "vep_cache.tar.gz")
        
        # Download into WORK_DIR
        run(f"aria2c -x 6 -s 6 -k 1M --check-certificate=false {ftp_url} -d {WORK_DIR} -o vep_cache.tar.gz")
        
        logger.info("Extracting cache with pigz...")
        run(f"tar -I pigz -xf {cache_tar} -C {VEP_DATA}")
        run(f"rm {cache_tar}")
    else:
        logger.info("Cache already present. Skipping download.")

    # 4. LOFTEE SETUP
    loftee_flags = ""
    if use_loftee:
        logger.info("Setting up LOFTEE resources from official repositories...")
        
        # 1. Plugin Scripts (Pull the entire grch38 branch from Konrad's repo)
        logger.info("Fetching LOFTEE grch38 scripts...")
        run(f"aria2c -x 5 -s 5 https://github.com/konradjk/loftee/archive/refs/heads/grch38.tar.gz -d {WORK_DIR} -o loftee_scripts.tar.gz")
        
        # --strip-components=1 unpacks the files directly into the Plugins folder without creating a subfolder
        run(f"tar -xzf {WORK_DIR}/loftee_scripts.tar.gz -C {PLUGINS_DIR} --strip-components=1")
        run(f"rm {WORK_DIR}/loftee_scripts.tar.gz")
        
        # 2. Data Files (From Konrad's official GRCh38 directory)
        broad_url = "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38"
        
        run(f"aria2c -x 5 -s 5 {broad_url}/loftee.sql.gz -d {VEP_DATA} -o loftee.sql.gz")
        run(f"pigz -d {VEP_DATA}/loftee.sql.gz") 
        
        # FIX: Force the file to be named exactly 'gerp.bw'
        run(f"aria2c -x 5 -s 5 {broad_url}/gerp_conservation_scores.homo_sapiens.GRCh38.bw -d {VEP_DATA} -o gerp.bw")
        
        # FIX: Force explicit names for the Ancestor files just to be 100% safe
        run(f"aria2c -x 5 -s 5 {broad_url}/human_ancestor.fa.gz -d {VEP_DATA} -o human_ancestor.fa.gz")
        run(f"aria2c -x 5 -s 5 {broad_url}/human_ancestor.fa.gz.fai -d {VEP_DATA} -o human_ancestor.fa.gz.fai")
        run(f"aria2c -x 5 -s 5 {broad_url}/human_ancestor.fa.gz.gzi -d {VEP_DATA} -o human_ancestor.fa.gz.gzi")

        loftee_flags = (
            f"--dir_plugins /opt/vep/.vep/Plugins "
            f"--plugin LoF,loftee_path:/opt/vep/.vep/Plugins,"
            f"human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,"
            f"conservation_file:/opt/vep/.vep/loftee.sql,"
            f"gerp_bigwig:/opt/vep/.vep/gerp.bw"
        )

    # 5. VEP EXECUTION
    CONTAINER_CACHE = "/opt/vep/.vep"
    docker_tag = f"ensemblorg/ensembl-vep:release_{vep_version}.0"
    logger.info(f"Pulling Docker image {docker_tag}...")
    run(f"docker pull {docker_tag}")

    logger.info("Opening permissions on vep_data so Docker can write output...")
    run(f"chmod -R 777 {VEP_DATA}")

    logger.info("Verify input is visible inside container via the cache mount...")
    run(f"docker run --rm -v {VEP_DATA}:{CONTAINER_CACHE} {docker_tag} ls -la {CONTAINER_CACHE}/input.vcf")

    vep_cmd = (
        f"docker run --rm "
        f"-v {VEP_DATA}:{CONTAINER_CACHE} "
        f"{docker_tag} vep "
        f"-i {CONTAINER_CACHE}/input.vcf -o {CONTAINER_CACHE}/output.tsv "
        f"--species homo_sapiens --assembly GRCh38 "
        f"--format vcf --tab --no_stats --cache --offline "
        f"--fork 12 "
        f"--af_gnomadg --af_gnomade "
        f"--total_length --no_escape "
        # f"--polyphen s "
        f"--canonical --protein --biotype "
        f"--dont_skip "
        f"--per_gene "
        f"--pick_order biotype,mane_select,canonical,appris,tsl,ccds,rank,length,ensembl,refseq "
        f"{loftee_flags}"
    )

    logger.info("Starting VEP annotation...")
    run(vep_cmd)

    # 6. FINAL UPLOAD
    local_output = os.path.join(VEP_DATA, "output.tsv")
    if not os.path.exists(local_output):
        raise FileNotFoundError(f"VEP finished but {local_output} was not created!")
        
    return {"chunk_tsv": dxpy.dxlink(dxpy.upload_local_file(local_output))}

@dxpy.entry_point("gather")
def gather(chunk_tsvs, use_loftee, master_parquet_link, genes_to_keep_file=None, biotypes_filter=None):
    
    logger.info("Downloading chunk results...")
    local_paths = []
    for idx, tsv_link in enumerate(chunk_tsvs):
        path = f"chunk_{idx}.tsv"
        dxpy.download_dxfile(tsv_link, path)
        local_paths.append(path)
        
    # --- DOWNLOAD MASTER PARQUET FOR THE JOIN ---
    logger.info("Downloading master variant metadata...")
    metadata_path = "variant_metadata.parquet"
    dxpy.download_dxfile(master_parquet_link, metadata_path)
        
    # --- Handle TXT Gene List ---
    gene_list = None
    if genes_to_keep_file:
        logger.info("Downloading custom gene list...")
        dxpy.download_dxfile(genes_to_keep_file, "custom_genes.txt")
        # Read the text file into a python list, stripping whitespace/newlines
        with open("custom_genes.txt", "r") as f:
            gene_list = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(gene_list)} genes for filtering.")
    
    # 1. Load all TSVs into a single LazyFrame
    logger.info("Concatenating VEP TSVs...")
    lazy_dfs = [
        pl.scan_csv(
            p,
            separator="\t",
            comment_prefix="##",
            null_values=["-", "."],
            infer_schema_length=0, # Highly recommended to prevent mid-stream crashes on VEP data
            ignore_errors=True,
        )
        for p in local_paths
    ]

    # Infer schema from first chunk, cast all gnomAD AF cols to Float64
    lf_colnames = lazy_dfs[0].collect_schema().names()
    schema_overrides = {
        col: pl.Float64
        for col in lf_colnames
        if col.startswith("gnomAD") and col.endswith("_AF")
    }

    lazy_dfs = [lf.with_columns(pl.col(col).cast(pl.Float64) for col in schema_overrides) for lf in lazy_dfs]
    full_lazy_df = pl.concat(lazy_dfs)
    
    # 2. APPLY POST-PROCESSING
    logger.info("Applying Post-Processing, Join, and LOFTEE transformations...")
    processed_lazy_df = post_process_vep(
        vep_lf=full_lazy_df,                  # Make sure this argument matches your updated post_process_vep
        metadata_parquet_path=metadata_path,  # Pass the path to the downloaded master parquet
        gene_list=gene_list,             
        biotypes_filter=biotypes_filter, 
        use_loftee=use_loftee
    )
    
    # 3. Stream the processed data to the final parquet file
    output_filename = "variants_vep_annotated.parquet"
    logger.info(f"Streaming final merged data to {output_filename}...")
    
    processed_lazy_df.sink_parquet(
        output_filename,
        engine='streaming'
    )
    
    logger.info("Gather step complete!")
    return {"vep_parquet": dxpy.dxlink(dxpy.upload_local_file(output_filename))}

@dxpy.entry_point("main")
def main(variants_parquet, chunk_size, vep_version, use_loftee, genes_to_keep_file=None, biotypes_filter=None):
    input_path = "input_main.parquet"
    dxpy.download_dxfile(variants_parquet, input_path)
    df = pl.read_parquet(input_path)
    logger.info(f"Total variants: {df.height}. Chunking by {chunk_size}...")

    subjob_outputs = []
    for i in range(0, df.height, chunk_size):
        chunk_path = f"chunk_{i}.parquet"
        df.slice(i, chunk_size).write_parquet(chunk_path)
        
        uploaded_file_obj = dxpy.upload_local_file(chunk_path)
        chunk_link = dxpy.dxlink(uploaded_file_obj)
        
        subjob = dxpy.new_dxjob(
            fn_name="process_chunk",
            fn_input={
                "chunk_file": chunk_link, 
                "vep_version": vep_version,
                "use_loftee": use_loftee
            },
            instance_type="mem2_ssd1_v2_x16"
        )
        subjob_outputs.append(subjob.get_output_ref("chunk_tsv"))
        
    gather_job = dxpy.new_dxjob(
        fn_name="gather", 
        fn_input={
            "chunk_tsvs": subjob_outputs, 
            "use_loftee": use_loftee,
            "master_parquet_link": variants_parquet,  # 2. Pass the input link straight down to gather
            "genes_to_keep_file": genes_to_keep_file, 
            "biotypes_filter": biotypes_filter        
        },
        instance_type="mem2_ssd1_v2_x16"
    )
    
    return {"vep_parquet": gather_job.get_output_ref("vep_parquet")}

dxpy.run()