#!/usr/bin/env python3

"""
TO BUILD: dx build applets/extract_phenotypes_and_covariates/ --destination <PROJECT-ID>:<PATH-TO-APPLET-DEST>
"""

import os
import subprocess
import dxpy
import polars as pl
from collections import defaultdict

PROJECT_ID = os.environ.get("DX_PROJECT_CONTEXT_ID")

def extract_phenotypes(file_link, dataset_id, batch_size=25, extra_fields=None):
    """
    Downloads a text file containing phenotype names, batches them to bypass 
    API limits, extracts the data via dx CLI, and returns a Polars DataFrame
    along with the list of extracted fields.
    """
    # Safely extract the file ID string from the DNAnexus link dictionary
    file_id = file_link if isinstance(file_link, str) else file_link["$dnanexus_link"]
    
    print(f"Downloading list file {file_id}...")
    local_list_name = f"list_{file_id}.txt"
    dxpy.download_dxfile(file_id, local_list_name)
    
    # Read the fields from the file
    with open(local_list_name, 'r') as f:
        fields_list = [line.strip() for line in f if line.strip()]
        
    if extra_fields:
        fields_list.extend(extra_fields)
        
    print(f"Found {len(fields_list)} fields. Preparing extraction...")

    # Clean, prefix, and deduplicate
    cleaned_fields = [f for f in fields_list if f.replace("participant.", "") != "eid"]
    prefixed_fields = [f"participant.{f}" if not f.startswith("participant.") else f for f in cleaned_fields]
    unique_fields = list(dict.fromkeys(prefixed_fields))

    full_dataset_path = f"{PROJECT_ID}:{dataset_id}"

    # Chunk into batches
    chunks = [unique_fields[i:i + batch_size] for i in range(0, len(unique_fields), batch_size)]
    chunk_dfs = []

    for i, chunk in enumerate(chunks):
        batch_fields = ["participant.eid"] + chunk
        fields_str = ",".join(batch_fields)
        out_file = f"chunk_{file_id}_{i}.tsv"
        
        print(f" -> Extracting Batch {i+1}/{len(chunks)} ({len(batch_fields)} columns)...")
        
        cmd = [
            "dx", "extract_dataset", full_dataset_path,
            "--fields", fields_str,
            "--output", out_file,
            "--delim", "\t" 
        ]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
            
            # Load into Polars and strip the "participant." prefix immediately
            df_chunk = pl.read_csv(out_file, separator="\t", null_values=["NA", ""])
            df_chunk = df_chunk.rename({col: col.replace("participant.", "") for col in df_chunk.columns})
            
            chunk_dfs.append(df_chunk)
            os.remove(out_file)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"FATAL: Batch {i+1} extraction failed. Check trait spelling! Error: {e}")

    # Merge all batches on 'eid'
    print("Merging batches into a single DataFrame...")
    final_df = chunk_dfs[0]
    for i, next_df in enumerate(chunk_dfs[1:]):
        assert next_df.height == chunk_dfs[0].height, f"Batch size mismatch at chunk {i}!"
        final_df = final_df.join(next_df, on="eid", how="inner")

    return final_df, cleaned_fields

def ukb_gen_related_with_data(data: pl.DataFrame, ukb_with_data: set, cutoff: float = 0.0884) -> pl.DataFrame:
    return data.filter(
        (pl.col("Kinship") > cutoff)
        & pl.col("ID1").is_in(ukb_with_data)
        & pl.col("ID2").is_in(ukb_with_data)
    )

def ukb_gen_samples_to_remove(data: pl.DataFrame, ukb_with_data: set, cutoff: float = 0.0884) -> list:
    data = ukb_gen_related_with_data(data, ukb_with_data=ukb_with_data, cutoff=cutoff)
    remove_samples: list = []

    while True:
        connections = (
            pl.concat([
                data.select(pl.col("ID1").alias("ID")),
                data.select(pl.col("ID2").alias("ID")),
            ])
            .group_by("ID")
            .agg(pl.len().alias("n"))
        )

        if connections.is_empty() or connections["n"].max() == 1:
            break

        most_connected_id = connections.sort(["n", "ID"], descending=[True, False]).row(0)[0]
        remove_samples.append(most_connected_id)

        data = data.filter(
            (pl.col("ID1") != most_connected_id) & (pl.col("ID2") != most_connected_id)
        )

    # remove_samples.extend(data["ID2"].to_list())
    remove_samples.extend(data["ID2"].unique().to_list())
    return remove_samples

def computed_unrelated_ids(df: pl.DataFrame, cutoff: float = 0.0884) -> pl.DataFrame:
    print("No unrelated file provided. Computing unrelated IDs from ukb_rel.dat...")
    
    # Download the kinship matrix directly using the DX CLI
    subprocess.run([
        "dx", "download",
        f"{PROJECT_ID}:/Bulk/Genotype Results/Genotype calls/ukb_rel.dat",
        "-o", "ukb_rel.dat"
    ], check=True)
    
    ukb_relatedness = pl.read_csv("ukb_rel.dat", separator=" ")
    current_samples_set = set(df["eid"].to_list())
    
    samples_to_remove = ukb_gen_samples_to_remove(
        ukb_relatedness, ukb_with_data=current_samples_set, cutoff=cutoff
    )
    
    unrelated_samples = list(current_samples_set - set(samples_to_remove))
    print(f"Kinship pruning complete. Removed {len(samples_to_remove)} related samples.")
    
    return pl.DataFrame({"eid": unrelated_samples})

def handle_multiple_measurements(df: pl.DataFrame, pheno_list: list) -> tuple[pl.DataFrame, list]:
    print("Checking for multiple measurements (within the same instance) to average...")
    base_to_cols = defaultdict(list)
    
    for col in pheno_list:
        if "_i" in col and "_a" in col and col in df.columns:
            base_instance = col.rsplit("_a", 1)[0]
            base_to_cols[base_instance].append(col)
            
    avg_exprs = []
    cols_to_drop = []
    new_phenos = []
    
    for base_instance, cols in base_to_cols.items():
        if len(cols) > 1:
            avg_exprs.append(pl.mean_horizontal(cols).alias(base_instance))
            cols_to_drop.extend(cols)
            new_phenos.append(base_instance)
            
    if avg_exprs:
        df = df.with_columns(avg_exprs).drop(cols_to_drop)
        pheno_list = [c for c in pheno_list if c not in cols_to_drop]
        pheno_list.extend(new_phenos)
        print(f" -> Averaged {len(cols_to_drop)} array columns into {len(new_phenos)} instance-specific columns.")
    else:
        print(" -> No multiple measurements found to average.")
        
    return df, pheno_list

def apply_statin_correction(df: pl.DataFrame) -> pl.DataFrame:
    STATIN_CODES = [
        "1140861958", # simvastatin
        "1140861970", # lipostat 10mg tablet
        "1140864592", # lescol 20mg capsule
        "1140881748", # zocor 10mg tablet
        "1140888594", # fluvastatin
        "1140888648", # pravastatin
        "1140910632", # eptastatin
        "1140910652", # synvinolin
        "1140910654", # velastatin
        "1141146138", # lipitor 10mg tablet
        "1141146234", # atorvastatin
        "1141192410", # rosuvastatin
        "1141192414", # crestor 10mg tablet
        "1141200040", # zocor heart-pro 10mg tablet
    ]
    
    print(f"Applying statin correction using {len(STATIN_CODES)} hardcoded codes...")
    med_cols = [c for c in df.columns if c.startswith("p20003")]
    if not med_cols:
        print("WARNING: No p20003 medication columns found! Skipping statin correction.")
        return df

    print("Evaluating medication arrays for statin usage...")
    df = df.with_columns(
        pl.any_horizontal([
            pl.col(c).cast(pl.String).is_in(STATIN_CODES) 
            for c in med_cols
        ]).fill_null(False).alias("uses_statin")
    )
    
    print(f" -> Identified {df['uses_statin'].sum()} participants taking statins.")

    chol_cols = [c for c in df.columns if c.startswith("p30690")]
    ldl_cols = [c for c in df.columns if c.startswith("p30780")]
    
    adjustments = []
    for col in chol_cols:
        adjustments.append(pl.when(pl.col("uses_statin")).then(pl.col(col) / 0.8).otherwise(pl.col(col)).alias(col))
    for col in ldl_cols:
        adjustments.append(pl.when(pl.col("uses_statin")).then(pl.col(col) / 0.7).otherwise(pl.col(col)).alias(col))
        
    if adjustments:
        df = df.with_columns(adjustments)
        print(f" -> Adjusted {len(chol_cols)} Cholesterol and {len(ldl_cols)} LDL columns.")

    return df

@dxpy.entry_point('main')
def main(dataset_id, pheno_list_file, covar_list_file, bim_file, subset_eur=False, subset_unrelated=False, statin_correction=False):
    
    # ---------------------------------------------------------
    # 1. DOWNLOAD INPUT FILES & EXTRACT DATA
    # ---------------------------------------------------------
    print("Extracting data from DNAnexus Apollo...")
    
    # If we are subsetting to EUR, we MUST extract the ancestry field (p30079)
    # We append it to covariates so we don't accidentally treat it as a trait later
    eur_field = ["p30079"] if subset_eur else []
    
    pheno_df, raw_pheno_list = extract_phenotypes(pheno_list_file, dataset_id)
    cov_df, raw_covar_list = extract_phenotypes(covar_list_file, dataset_id, extra_fields=eur_field)
    
    print("Merging phenotype and covariate dataframes...")
    df = pheno_df.join(cov_df, on="eid", how="inner")
    print(f"Total merged dataset: {df.height} participants, {len(df.columns)} columns.")

    print("Downloading bim file...")
    # Safely extract the file ID string from the DNAnexus link dictionary
    bim_id = bim_file if isinstance(bim_file, str) else bim_file["$dnanexus_link"]
    dxpy.download_dxfile(bim_id, "input.bim")

    # ---------------------------------------------------------
    # 2. RUN QC AND FILTERING 
    # ---------------------------------------------------------
    
    # A. Average Multiple Measurements
    df, final_pheno_list = handle_multiple_measurements(df, raw_pheno_list)

    # B. Subset to EUR Ancestry
    if subset_eur:
        print("Subsetting to EUR ancestry...")
        # Find the ancestry column dynamically (handles p30079 or p30079_i0)
        eur_col = [c for c in df.columns if c.startswith("p30079")]
        if eur_col:
            df = df.filter(pl.col(eur_col[0]) == 1)             # EUR code is 1 in field 30079
            print(f"Retained {df.height} EUR participants.")
        else:
            print("WARNING: Ancestry field 30079 not found. Skipping EUR subset.")

    # C. Subset to Unrelated
    if subset_unrelated:
        print("Subsetting to unrelated individuals...")
        unrelated_df = computed_unrelated_ids(df)
        df = df.join(unrelated_df, on="eid", how="semi")
        print(f"Retained {df.height} unrelated participants.")

    # D. Statin Correction
    if statin_correction:
        df = apply_statin_correction(df)

    # ---------------------------------------------------------
    # 3. FORMAT FOR REGENIE
    # ---------------------------------------------------------
    print("Formatting FID and IID columns for REGENIE...")
    df = df.with_columns([
        pl.col("eid").alias("FID"),
        pl.col("eid").alias("IID")
    ])

    # Safely pull the requested columns that survived the QC pipeline
    final_covar_cols = ["FID", "IID"] + [c for c in raw_covar_list if c in df.columns]
    final_pheno_cols = ["FID", "IID"] + [c for c in final_pheno_list if c in df.columns]

    covar_out = "covariates.txt"
    pheno_out = "phenotypes.txt"

    df.select(final_covar_cols).write_csv(covar_out, separator=" ", null_value="NA")
    df.select(final_pheno_cols).write_csv(pheno_out, separator=" ", null_value="NA")
    print(f"Wrote covariates ({len(final_covar_cols)-2} cols) and phenotypes ({len(final_pheno_cols)-2} cols).")

    # ---------------------------------------------------------
    # 4. PROCESS BIM FILE TO GENERATE SNPLIST
    # ---------------------------------------------------------
    print("Extracting SNP list from the .bim file...")
    # .bim files are 6 columns with no header. The variant ID is the 2nd column.
    bim_df = pl.read_csv("input.bim", separator="\t", has_header=False, 
                         new_columns=["CHR", "SNP", "CM", "BP", "A1", "A2"], infer_schema_length=0)
    
    snplist_out = "snplist.snplist"
    # Write only the SNP column, no header, for REGENIE
    bim_df.select("SNP").write_csv(snplist_out, has_header=False)
    print(f"Extracted {bim_df.height} variants into snplist.")

    # ---------------------------------------------------------
    # 5. UPLOAD OUTPUTS BACK TO DNANEXUS
    # ---------------------------------------------------------
    print("Uploading output files back to DNAnexus...")
    return {
        "covariates": dxpy.dxlink(dxpy.upload_local_file(covar_out)),
        "phenotypes": dxpy.dxlink(dxpy.upload_local_file(pheno_out)),
        "snplist": dxpy.dxlink(dxpy.upload_local_file(snplist_out))
    }

dxpy.run()