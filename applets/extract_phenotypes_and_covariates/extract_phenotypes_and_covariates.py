#!/usr/bin/env python3

"""
TO BUILD: dx build applets/extract_phenotypes_and_covariates/ --destination <PROJECT-ID>:<PATH-TO-APPLET-DEST>
"""

import os
import subprocess
import dxpy
import polars as pl
from collections import defaultdict
from scipy.special import ndtri
from sklearn.linear_model import LinearRegression
import re

PROJECT_ID = os.environ.get("DX_PROJECT_CONTEXT_ID")

def extract_phenotypes(fields_list, dataset_id, batch_size=25):
    """
    Downloads phenotype fields speciffied in the fields_list, batches them to bypass 
    API limits, extracts the data via dx CLI, and returns a Polars DataFrame
    along with the list of extracted fields.
    """
    print(f"Preparing extraction of {len(fields_list)} fields...")

    # Clean, prefix, and deduplicate
    cleaned_fields = [f for f in fields_list if f.replace("participant.", "") != "eid"]
    prefixed_fields = [f"participant.{f}" if not f.startswith("participant.") else f for f in cleaned_fields]
    unique_fields = list(dict.fromkeys(prefixed_fields))

    dataset_str = dataset_id if isinstance(dataset_id, str) else dataset_id["$dnanexus_link"]
    full_dataset_path = f"{PROJECT_ID}:{dataset_str}"

    # Chunk into batches
    chunks = [unique_fields[i:i + batch_size] for i in range(0, len(unique_fields), batch_size)]
    chunk_dfs = []

    for i, chunk in enumerate(chunks):
        batch_fields = ["participant.eid"] + chunk
        fields_str = ",".join(batch_fields)
        out_file = f"pheno_chunk_{i}.tsv"
        
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
    
    # Extract field labels
    field_labels_out = 'field_labels.tsv'
    cmd = ["dx", "extract_dataset", full_dataset_path, "--list-fields",]
            
    with open(field_labels_out, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)
    
    field_labels = (
        pl.read_csv(field_labels_out, separator = "\t", has_header=False)
        .rename({'column_1': 'field_id', 'column_2': 'field_label'})
        .filter(pl.col('field_id').is_in(unique_fields))
    )

    # Create a renaming dict
    field_renaming = {}
    for row in field_labels.iter_rows():
        field = row[0].replace("participant.", "")
        label = row[1].split('|')[0].strip().lower()

        # Remove anything that isn't a letter, number, or space
        label = re.sub(r'[^a-z0-9\s]', '', label)

        # Replace one or more spaces with a single underscore
        label = re.sub(r'\s+', '_', label)

        # Preserve only the array suffix (_aN) so array fields that share
        # a single label (p22009_a1..a20, p20003_a0..a47, ...) don't
        # collapse into duplicate column names
        m = re.search(r'(_a\d+)$', field)
        suffix = m.group(1) if m else ""
        field_renaming[field] = f"{label}{suffix}"

    # Remove "participant." from field names
    unique_fields = [i.replace('participant.', "") for i in unique_fields]
    
    return final_df, unique_fields, field_renaming

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

def handle_multiple_measurements(df: pl.DataFrame, pheno_list: list, field_labelling: dict):
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
            # Inherit the label from the first array col, but drop the _aN
            # suffix since the averaged column is no longer array-specific
            inherited_label = re.sub(r'_a\d+$', '', field_labelling[cols[0]])
            field_labelling[base_instance] = inherited_label
            # remove old field labels
            for c in cols:
                del field_labelling[c]
            
    if avg_exprs:
        df = df.with_columns(avg_exprs).drop(cols_to_drop)
        pheno_list = [c for c in pheno_list if c not in cols_to_drop]
        pheno_list.extend(new_phenos)
        print(f" -> Averaged {len(cols_to_drop)} array columns into {len(new_phenos)} instance-specific columns.")
    else:
        print(" -> No multiple measurements found to average.")
        
    return df, pheno_list, field_labelling

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

def apply_quantile_transform(df: pl.DataFrame, phenotype_list: list) -> pl.DataFrame:
    
    c = 3/8  # Blom's constant for inverse normal transformation (prevents infinite values at the tails)
    
    # Create a long phenotypes files
    long_phenos_int = (
        df
        .select(['eid'] + phenotype_list)
        .unpivot(
            index='eid',
            on=phenotype_list,
            variable_name='phenotype',
            value_name='pheno_value'
        )
        .lazy()  # Use Lazy mode for better memory/query optimization
        .with_columns(
            # Calculate rank and group size using native Rust engine
            r = pl.col("pheno_value").rank().over("phenotype"),
            n = pl.len().over("phenotype")
        )
        .with_columns(
            # Calculate the INT value calling ndtri ONCE on the whole column
            pheno_value_int = ((pl.col("r") - c) / (pl.col("n") - 2*c + 1)).map_batches(ndtri)
        )
        .drop(["r", "n"]) # Clean up temporary columns
        .collect()
    )

    # Transform long to wide
    wide_phenos = (
        long_phenos_int
        .drop('pheno_value')
        .pivot(
        on="phenotype",           
        index="eid",              
        values="pheno_value_int")
    )
    
    # Add QT-corrected values to df
    df = (df
        .drop(phenotype_list)
        .join(wide_phenos, on = 'eid', how = 'inner')
    )

    return df

def process_regenie_ouputs(PRS_out_filename: str, PRS_out_folder: str) -> pl.DataFrame:

    # Read a mapping between between phenotype name and a filename
    PRS_phenos_to_files = (
        pl.read_csv(f"{PRS_out_filename}_prs.list", has_header = False, separator = " ")
        .rename({'column_1': 'phenotype', 'column_2': 'filename'})
    )

    # Read all PRS files
    prs_outputs = []
    for row in PRS_phenos_to_files.iter_rows():
        phenotype = row[0]
        filename = row[1]
        prs_one_pheno = (
            pl.read_csv(f"{PRS_out_folder}{filename}", separator = " ", null_values = "NA")
            .transpose(
                include_header=True, 
                header_name="eid",     
                column_names=['prs_value']  
            )
            .filter(pl.col("eid") != "FID_IID")
            .filter(pl.col("eid") != "")
            .with_columns(
                eid = pl.col("eid").str.split("_").list.get(0).cast(pl.Int64),
                prs_value = pl.col('prs_value').cast(pl.Float64)
            )
            .rename({'prs_value': f"{phenotype}_prs"})
        )
        prs_outputs.append(prs_one_pheno)
    
    # Collect all PRS files into a single dataframe
    PRS_df = prs_outputs[0]
    for i, next_df in enumerate(prs_outputs[1:]):
        PRS_df = PRS_df.join(next_df, on="eid", how="inner")

    return PRS_df

def correct_phenotypes(df: pl.DataFrame, phenotype_list: list, covariates_list: list) -> pl.DataFrame:
    
    # Correct phenos one by one
    correted_dfs = []
    for pheno_col in phenotype_list:

        # Add corresponding PRS to feature columns
        feature_cols = covariates_list + [f'{pheno_col}_prs']

        # Drop null values to residualize
        df_per_pheno = (
            df
            .select(['eid', pheno_col] + feature_cols)
            .filter(pl.col(pheno_col).is_not_null())
            )

        # Convert to NumPy for sklearn 
        X = df_per_pheno.select(feature_cols).to_numpy()
        y = df_per_pheno.select(pheno_col).to_numpy()

        # Fit the model and get residuals
        model = LinearRegression()
        model.fit(X, y)
        residuals = y - model.predict(X)
        
        # Put residuals into a Dataframe
        df_per_pheno = (
            df_per_pheno
            .select(['eid'])
            .with_columns(pl.Series(residuals.flatten()).cast(pl.Float64).alias(f"{pheno_col}"))
        )

        # Add null values back
        df_per_pheno = (
            df
            .select(['eid'])
            .join(df_per_pheno, on = 'eid', how = 'left', coalesce = True)
        )
        correted_dfs.append(df_per_pheno) 

    # Collect all residuals into a single dataframe
    corrected_df = correted_dfs[0]
    for i, next_df in enumerate(correted_dfs[1:]):
        corrected_df = corrected_df.join(next_df, on="eid", how="inner")    
    
    return corrected_df
        
@dxpy.entry_point('main')
def main(dataset_id, 
        pheno_list_file, 
        bed_file=None, 
        bim_file=None, 
        fam_file=None, 
        additional_covars_file=None, 
        subset_eur=False, 
        subset_unrelated=False, 
        statin_correction=False):
    
    # ---------------------------------------------------------
    # 1. DOWNLOAD INPUT FILES & EXTRACT DATA
    # ---------------------------------------------------------
    print("Extracting data from DNAnexus Apollo...")
    
    # Define standard REGENIE covariates: Age (p21003_i0), Sex (p22001), and PCs 1-20
    standard_covars = ["p21003_i0", "p22001"] + [f"p22009_a{i}" for i in range(1, 21)]

    eur_field = ["p30079"] if subset_eur else []
    medication_fields = [f'p20003_i0_a{i}' for i in range(0, 48)] if statin_correction else []
    # Append any user-provided extra covariates
    if additional_covars_file:
        additional_covars_file_id = additional_covars_file if isinstance(additional_covars_file, str) else additional_covars_file["$dnanexus_link"]
        print(f"Downloading list file {additional_covars_file_id}...")
        local_list_name = f"list_{additional_covars_file_id}.txt"
        dxpy.download_dxfile(additional_covars_file_id, local_list_name)
        
        # Read the fields from the file
        with open(local_list_name, 'r') as f:
            additional_covars_fields = [line.strip() for line in f if line.strip()]
    else:
        additional_covars_fields = []

    covar_fields_list = standard_covars + eur_field + medication_fields + additional_covars_fields

    cov_df, raw_cov_list, covariates_field_labelling = extract_phenotypes(covar_fields_list, dataset_id)
    # remove ancestry and medication from the covariates list
    final_cov_list = [i for i in raw_cov_list if i not in medication_fields]
    final_cov_list = [i for i in final_cov_list if i not in eur_field]

    # Safely extract the file ID string from the DNAnexus link dictionary
    file_id = pheno_list_file if isinstance(pheno_list_file, str) else pheno_list_file["$dnanexus_link"]
    print(f"Downloading list file {file_id}...")
    local_list_name = f"list_{file_id}.txt"
    dxpy.download_dxfile(file_id, local_list_name)
    
    # Read the fields from the file
    with open(local_list_name, 'r') as f:
        pheno_fields_list = [line.strip() for line in f if line.strip()]

    pheno_df, raw_pheno_list, phenotypes_field_labelling = extract_phenotypes(pheno_fields_list, dataset_id)
    
    print("Merging phenotype and covariate dataframes...")
    df = pheno_df.join(cov_df, on="eid", how="inner")
    print(f"Total merged dataset: {df.height} participants, {len(df.columns)} columns.")

    corrected_phenos_pq_out = "corrected_phenos.parquet"
    phenos_pq_out = "phenos.parquet"
    covariates_pq_out = "covariates.parquet"
    prs_pq_out = "prs.parquet"
    
    # ---------------------------------------------------------
    # 2. RUN QC AND FILTERING 
    # ---------------------------------------------------------
    
    # A. Average Multiple Measurements
    df, final_pheno_list, phenotypes_field_labelling = handle_multiple_measurements(df, raw_pheno_list, phenotypes_field_labelling)

    # B. Subset to EUR Ancestry
    if subset_eur:
        print("Subsetting to EUR ancestry...")
        # Find the ancestry column dynamically (handles p30079 or p30079_i0)
        eur_col = [c for c in df.columns if c.startswith("p30079")]
        if eur_col:
            df = df.filter(pl.col(eur_col[0]) == 5)             # EUR code is 5 in field 30079
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

    # E. Apply quantile transformation
    df = apply_quantile_transform(df, final_pheno_list)

    # ---------------------------------------------------------
    # 3. FORMAT FOR REGENIE
    # ---------------------------------------------------------
    print("Formatting FID and IID columns for REGENIE...")
    df = df.with_columns([
        pl.col("eid").alias("FID"),
        pl.col("eid").alias("IID")
    ])

    # Safely pull the requested columns that survived the QC pipeline ()
    final_covar_cols = ["FID", "IID"] + [c for c in final_cov_list if c in df.columns]
    final_pheno_cols = ["FID", "IID"] + [c for c in final_pheno_list if c in df.columns]

    covar_out = "covariates.txt"
    pheno_out = "phenotypes.txt"
    
    # Save phenos and covariates for regenie
    df.select(final_covar_cols).write_csv(covar_out, separator=" ", null_value="NA")
    df.select(final_pheno_cols).write_csv(pheno_out, separator=" ", null_value="NA")

    # Return to normal identifier format
    # df = df.drop('FID').rename({'IID': 'eid'})
    
    # Save parquet covariates
    (
        df
        .select(['eid'] + [c for c in final_cov_list if c in df.columns])
        .rename(covariates_field_labelling, strict=False)
        .write_parquet(covariates_pq_out)
    )

    # Save parquet phenotypes
    (
        df
        .select(['eid'] + [c for c in final_pheno_list if c in df.columns])
        .rename(phenotypes_field_labelling, strict=False)
        .write_parquet(phenos_pq_out)
    )

    print(f"Wrote covariates ({len(final_covar_cols)-2} cols) and phenotypes ({len(final_pheno_cols)-2} cols).")

    # ---------------------------------------------------------
    # 4. RUN REGENIE IF USER PROVIDED .bed, .bim and .fam FILES 
    # ---------------------------------------------------------
    if bed_file and bim_file and fam_file:
        genotype_calls_out = "input"
        PRS_out_folder = "./regenie_out/"
        PRS_out_filename = f"{PRS_out_folder}file"
        os.makedirs(PRS_out_folder, exist_ok=True)

        print("Downloading genotype call files...")
        # Safely extract the file ID string from the DNAnexus link dictionary
        bed_id = bed_file if isinstance(bed_file, str) else bed_file["$dnanexus_link"]
        dxpy.download_dxfile(bed_id, f"{genotype_calls_out}.bed")

        bim_id = bim_file if isinstance(bim_file, str) else bim_file["$dnanexus_link"]
        dxpy.download_dxfile(bim_id, f"{genotype_calls_out}.bim")

        fam_id = fam_file if isinstance(fam_file, str) else fam_file["$dnanexus_link"]
        dxpy.download_dxfile(fam_id, f"{genotype_calls_out}.fam")

        print("Installing conda environment...")

        subprocess.run("wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh", shell=True)
        subprocess.run("bash miniconda.sh -b -p $HOME/miniconda", shell=True)
        # Add conda to the current path
        conda_bin = os.path.expanduser("~/miniconda/bin/conda")

        print("Creating REGENIE environment (this takes 3-5 minutes)...")
        subprocess.run([
            conda_bin, "create", "-y", "-n", "regenie_env", 
            "--override-channels", 
            "-c", "conda-forge", 
            "-c", "bioconda", 
            "regenie"
        ], check=True)

        print("Starting REGENIE Step 1...")
        regenie_args = [
            "regenie",
            "--step", "1",
            # Genotypes (to .bed without ".bed")
            "--bed", genotype_calls_out,

            # Phenotypes & Covariates
            "--phenoFile", pheno_out,
            "--covarFile", covar_out,

            # Model Parameters
            "--bsize", "1000",      # normally we use this parameter
            "--qt",                 # indicate that phenos are quantitative
            "--lowmem",             # Low-memory to not crash the session

            # PRS Mode & Output
            "--print-prs",          # Force output of .prs files (Whole-Genome Predictions)
            "--out", PRS_out_filename
        ]

        full_cmd = [conda_bin, "run", "-n", "regenie_env"] + regenie_args

        try:
            result = subprocess.run(full_cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print("--- REGENIE FAILED ---")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr) 
            raise
        
        print("Processing regenie outputs...")
        PRS_df = process_regenie_ouputs(PRS_out_filename, PRS_out_folder)
        df = PRS_df.join(df, on = 'eid', how = 'inner')

        # Save PRS: strip _prs suffix from PRS columns so phenotypes_field_labelling applies
        (
            PRS_df
            .rename({f"{i}_prs": i for i in final_pheno_list}, strict=False)
            .rename(phenotypes_field_labelling, strict=False)
            .write_parquet(prs_pq_out)
        )
        
        
        # ---------------------------------------------------------
        # 5. CORRECT FOR COVARIATES AND PRS 
        # ---------------------------------------------------------
        corrected_phenotypes = correct_phenotypes(df, final_pheno_list, final_cov_list)
        corrected_phenotypes.rename(phenotypes_field_labelling, strict=False).write_parquet(corrected_phenos_pq_out)

    else:
        print("No genotype call files provided, skipping PRS calculation and phenotype correction steps.")
        corrected_phenotypes = None
        

    # ---------------------------------------------------------
    # 6. UPLOAD OUTPUTS BACK TO DNANEXUS
    # ---------------------------------------------------------

    print("Uploading output files back to DNAnexus...")
    # RETURN: covariates (remove medication), qt-phenos, PRS and adjusted phenos
    return {
        "covariates": dxpy.dxlink(dxpy.upload_local_file(covariates_pq_out)),
        "phenotypes": dxpy.dxlink(dxpy.upload_local_file(phenos_pq_out)),
        "PRS": dxpy.dxlink(dxpy.upload_local_file(prs_pq_out)) if bed_file and bim_file and fam_file else None,
        "corrected_phenotypes": dxpy.dxlink(dxpy.upload_local_file(corrected_phenos_pq_out)) if bed_file and bim_file and fam_file else None,
    }

dxpy.run()