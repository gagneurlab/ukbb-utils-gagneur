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
import pyranges as pr

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Single source of truth for paths
WORK_DIR = "/home/dnanexus"
VEP_DATA = os.path.join(WORK_DIR, "vep_data")

# Download URLs
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/"
    "gencode.v39.primary_assembly.annotation.gtf.gz"
)
REFERENCE_FASTA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/"
    "GRCh38.primary_assembly.genome.fa.gz"
)

def run(cmd):
    logger.info(f"$ {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def aria2c_download(url: str, outpath: str) -> bool:
    """Download a file with aria2c. Returns True on success, False on failure."""
    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    try:
        run(
            f"aria2c -x 16 -s 16 --allow-overwrite=true "
            f"'{url}' -d '{os.path.dirname(outpath)}' -o '{os.path.basename(outpath)}'"
        )
        return os.path.exists(outpath)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Download failed for {url}: {e}")
        return False

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

def _convert_to_int_and_get_max(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        if value is None:
            return None
        parts = str(value).split("-")
        ints = [int(v) for v in parts if v.isdigit()]
        return max(ints) if ints else None


def _next_inframe_start_codon_distance(seq, search_init=0):
    pos = seq.find("ATG", search_init)
    if pos == -1:
        return -1
    rel = pos - search_init
    if rel % 3 == 0:
        return pos
    next_init = pos + (3 - rel % 3)
    if next_init >= len(seq):
        return -1
    return _next_inframe_start_codon_distance(seq, next_init)


def _compute_next_in_frame(annos: pl.LazyFrame, fasta_path: str) -> pl.LazyFrame:
    """Compute next in-frame ATG distance for start_lost variants (requires FASTA)."""
    import pandas as pd

    schema_names = set(annos.collect_schema().names())
    req_cols = {"pos", "chrom", "region", "id", "cds_position", "codons", "strand", "allele"}
    if not req_cols.issubset(schema_names):
        logger.warning("  Missing columns for next_in_frame; skipping")
        return annos

    vep_start_lost = (
        annos.filter(pl.col("consequence_start_lost") == 1)
        .select(list(req_cols))
        .collect()
        .to_pandas()
    )
    snvs = vep_start_lost[vep_start_lost["allele"].str.len() == 1].copy()
    snvs[["variant_pos_cds", "cds_length"]] = snvs["cds_position"].str.split("/", expand=True)
    snvs = snvs[snvs["variant_pos_cds"].str.len() == 1]
    snvs["variant_pos_cds"] = snvs["variant_pos_cds"].astype(int)
    snvs["cds_length"] = snvs["cds_length"].astype(int)
    snvs["Chromosome"] = snvs["chrom"].str.replace(r"^(?!chr)", "chr", regex=True)

    results = []
    for strand_val in ["1", "-1"]:
        sub = snvs[snvs["strand"] == strand_val].copy()
        if sub.empty:
            continue
        sub["Start"] = (
            sub["pos"] - sub["variant_pos_cds"]
            if strand_val == "1"
            else sub["pos"] + sub["variant_pos_cds"] - 1 - sub["cds_length"]
        )
        sub["End"] = sub["Start"] + sub["cds_length"] + 100
        sub["Strand"] = "+" if strand_val == "1" else "-"
        sub_pr = pr.PyRanges(sub)
        sub_pr.seq = pr.get_sequence(sub_pr, path=fasta_path)
        sub = sub_pr.df.copy()
        sub["next_in_frame"] = sub.apply(
            lambda x: _next_inframe_start_codon_distance(x["seq"], search_init=3), axis=1
        )
        sub["next_in_frame_relative"] = sub["next_in_frame"] / sub["cds_length"]
        sub.loc[sub["next_in_frame_relative"] < 0, "next_in_frame_relative"] = 1
        results.append(sub[["region", "id", "next_in_frame_relative"]])

    if results:
        nif_pl = pl.from_pandas(pd.concat(results))
        annos = annos.join(nif_pl.lazy(), on=["id", "region"], how="left", validate="1:1")
    return annos


def _parse_vep_extra_col(lf: pl.LazyFrame, use_polyphen: bool = False) -> pl.LazyFrame:
    """
    Extracts LOFTEE (and optionally PolyPhen2) fields from the VEP 'extra' column.
    """
    logger.info("Extracting LOFTEE and standard fields from 'extra' column...")

    if "extra" not in lf.collect_schema().names():
        return lf

    extra_exprs = [
        pl.col("extra").str.extract(r"(?:^|;)LoF=([^;]*)", 1).alias("lof"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_filter=([^;]*)", 1).alias("lof_filter"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_flags=([^;]*)", 1).alias("lof_flags"),
        pl.col("extra").str.extract(r"(?:^|;)LoF_info=([^;]*)", 1).alias("lof_info"),
        pl.col("extra").str.extract(r"(?:^|;)IMPACT=([^;]*)", 1).alias("impact"),
    ]

    if use_polyphen:
        # PolyPhen2 format: PolyPhen=possibly_damaging(0.876)
        extra_exprs += [
            pl.col("extra").str.extract(r"(?:^|;)PolyPhen=([^(;]+)", 1).alias("polyphen_prediction"),
            pl.col("extra").str.extract(r"(?:^|;)PolyPhen=[^(]*\(([^)]*)\)", 1)
                .cast(pl.Float32).alias("polyphen_score"),
        ]

    return lf.with_columns(extra_exprs).drop("extra")

def add_vep_structural_features(
    annos: pl.LazyFrame,
    gtf_path: str,
    ref_fasta_path: str = None,
) -> pl.LazyFrame:
    """
    Add VEP-derived structural features: relative_cds_position, variant_length/indel flags,
    dist_to_tss, gene_length, gene_name (from GTF), next_in_frame_relative (if FASTA provided).
    """
    schema = set(annos.collect_schema().names())

    # ── Relative CDS position ─────────────────────────────────────────────
    logger.info("  Adding relative CDS position")
    if "cds_position" in schema:
        try:
            cds_parsed = (
                annos.with_columns(
                    pl.col("cds_position").str.split("/").alias("_cds_parts")
                )
                .with_columns(
                    _length=pl.col("_cds_parts").list.get(1),
                    _prot_pos_str=pl.col("_cds_parts").list.get(0),
                )
                .with_columns(
                    _prot_pos_int=pl.col("_prot_pos_str").map_elements(
                        _convert_to_int_and_get_max, return_dtype=pl.Int64
                    )
                )
                .filter(pl.col("_prot_pos_int").is_not_null())
                .with_columns(pl.col("_length").cast(pl.Int64, strict=False))
                .with_columns(
                    relative_cds_position=(
                        pl.col("_prot_pos_int") / pl.col("_length")
                    ).round(2)
                )
                .select(["id", "region", "relative_cds_position"])
            )
            annos = annos.join(cds_parsed, on=["id", "region"], how="left")
        except Exception as e:
            logger.warning(f"    relative_cds_position failed: {e}")

    # ── Variant length and indel flags ────────────────────────────────────
    logger.info("  Adding variant length and indel flags")
    annos = annos.with_columns(
        pl.max_horizontal(
            [pl.col("ref").str.len_chars(), pl.col("alt").str.len_chars()]
        ).alias("variant_length")
    ).with_columns(
        is_indel=(pl.col("variant_length") > 1).cast(pl.Int8),
        is_insertion=(
            pl.col("ref").str.len_chars() < pl.col("alt").str.len_chars()
        ).cast(pl.Int8),
        is_deletion=(
            pl.col("ref").str.len_chars() > pl.col("alt").str.len_chars()
        ).cast(pl.Int8),
    )

    # ── Distance to TSS, gene_length, gene_name from Gencode GTF ──────────
    logger.info("  Adding dist_to_tss, gene_length, gene_name from GTF")
    try:
        gencode_pr = pr.read_gtf(gtf_path, as_df=True)
        gencode_pl = pl.from_pandas(gencode_pr).filter(
            pl.col("gene_type") == "protein_coding"
        )
        gencode_genes = gencode_pl.filter(pl.col("Feature") == "gene").with_columns(
            region=pl.col("gene_id").str.split(".").list.first()
        )
        tss_df = gencode_genes.with_columns(
            gene_length=pl.col("End") - pl.col("Start") + 1,
            tss=pl.when(pl.col("Strand") == "+")
            .then(pl.col("Start"))
            .otherwise(pl.col("End")),
        ).select(["region", "gene_name", "gene_length", "tss", "Strand"])

        anno_tss = (
            annos.select(["id", "pos", "region"])
            .join(tss_df.lazy(), on="region", how="left")
            .with_columns(
                dist_to_tss=pl.when(pl.col("Strand") == "+")
                .then(pl.col("pos") - pl.col("tss"))
                .otherwise(pl.col("tss") - pl.col("pos"))
            )
            .select(["id", "region", "dist_to_tss", "gene_length", "gene_name"])
        )
        annos = annos.join(anno_tss, on=["id", "region"], how="left")
        logger.info("    dist_to_tss computed successfully")
    except Exception as e:
        logger.warning(f"    dist_to_tss failed: {e}")

    # ── next_in_frame_relative (start_lost, requires FASTA) ───────────────
    if ref_fasta_path and "consequence_start_lost" in annos.collect_schema().names():
        logger.info("  Computing next_in_frame_relative for start_lost variants")
        try:
            annos = _compute_next_in_frame(annos, ref_fasta_path)
        except Exception as e:
            logger.warning(f"    next_in_frame_relative failed: {e}")

    return annos


def post_process_vep(vep_lf, metadata_parquet_path, gene_list=None, biotypes_filter=None, use_loftee=True, use_polyphen=False):
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
    lf = _parse_vep_extra_col(lf, use_polyphen=use_polyphen)

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
def process_chunk(chunk_file, vep_version, use_loftee, use_polyphen=False):
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

    polyphen_flag = "--polyphen s " if use_polyphen else ""

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
        f"{polyphen_flag}"
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
def gather(chunk_tsvs, use_loftee, master_parquet_link, genes_to_keep_file=None, biotypes_filter=None, use_polyphen=False):

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

    # --- DOWNLOAD GENCODE GTF AND REFERENCE FASTA ---
    logger.info("Downloading Gencode v39 GTF...")
    gtf_path = os.path.join(WORK_DIR, "gencode.v39.primary_assembly.annotation.gtf.gz")
    ok_gtf = aria2c_download(GENCODE_GTF_URL, gtf_path)
    if not ok_gtf:
        logger.warning("GTF download failed; structural features will be skipped")
        gtf_path = None

    logger.info("Downloading Gencode v39 reference FASTA...")
    fasta_path = os.path.join(WORK_DIR, "GRCh38.primary_assembly.genome.fa.gz")
    ok_fasta = aria2c_download(REFERENCE_FASTA_URL, fasta_path)
    if not ok_fasta:
        logger.warning("FASTA download failed; next_in_frame features will be skipped")
        fasta_path = None

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
        vep_lf=full_lazy_df,
        metadata_parquet_path=metadata_path,
        gene_list=gene_list,
        biotypes_filter=biotypes_filter,
        use_loftee=use_loftee,
        use_polyphen=use_polyphen,
    )

    # 3. ADD VEP STRUCTURAL FEATURES
    if gtf_path:
        logger.info("Adding VEP structural features (CDS position, indel flags, TSS distance, gene annotations)...")
        processed_lazy_df = add_vep_structural_features(
            annos=processed_lazy_df,
            gtf_path=gtf_path,
            ref_fasta_path=fasta_path if fasta_path else None,
        )
    else:
        logger.warning("Skipping structural features (GTF unavailable)")
    
    # 4. Stream the processed data to the final parquet file
    output_filename = "variants_vep_annotated.parquet"
    logger.info(f"Streaming final merged data to {output_filename}...")

    processed_lazy_df.sink_parquet(
        output_filename,
        engine='streaming'
    )

    logger.info("Gather step complete!")
    return {"vep_parquet": dxpy.dxlink(dxpy.upload_local_file(output_filename))}

@dxpy.entry_point("main")
def main(variants_parquet, chunk_size, vep_version, use_loftee, use_polyphen=False, genes_to_keep_file=None, biotypes_filter=None):
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
                "use_loftee": use_loftee,
                "use_polyphen": use_polyphen,
            },
            instance_type="mem2_ssd1_v2_x16"
        )
        subjob_outputs.append(subjob.get_output_ref("chunk_tsv"))

    gather_job = dxpy.new_dxjob(
        fn_name="gather",
        fn_input={
            "chunk_tsvs": subjob_outputs,
            "use_loftee": use_loftee,
            "use_polyphen": use_polyphen,
            "master_parquet_link": variants_parquet,
            "genes_to_keep_file": genes_to_keep_file,
            "biotypes_filter": biotypes_filter,
        },
        instance_type="mem2_ssd1_v2_x16"
    )
    
    return {"vep_parquet": gather_job.get_output_ref("vep_parquet")}

dxpy.run()