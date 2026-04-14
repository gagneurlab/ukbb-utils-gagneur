#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polars",
# ]
# ///

import argparse
from pathlib import Path

import polars as pl


def main():
    parser = argparse.ArgumentParser(description="Prepare regenie Olink input files.")
    parser.add_argument("--covariates-parquet", required=True, type=Path)
    parser.add_argument("--levels-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- proteomics levels ---
    proteomics_levels = pl.read_csv(args.levels_csv)
    genes = proteomics_levels["gene"].to_list()
    proteomics_levels = (
        proteomics_levels.select(pl.exclude("gene"))
        .transpose(header_name="sample", include_header=True, column_names=genes)
        .rename({"sample": "FID"})
        .with_columns(FID=pl.col("FID").cast(pl.Int32))
        .with_columns(IID=pl.col("FID"))
        .select(["FID", "IID"] + genes)
    )

    # --- covariates ---
    proteomics_covariates = pl.read_parquet(args.covariates_parquet)
    covariate_columns = [c for c in proteomics_covariates.columns if c != "sample"]
    proteomics_covariates = (
        proteomics_covariates.rename({"sample": "FID"})
        .with_columns(FID=pl.col("FID").cast(pl.Int32))
        .with_columns(IID=pl.col("FID"))
        .select(["FID", "IID"] + covariate_columns)
        .drop_nulls()
        .sort("FID")
    )
    proteomics_covariates.columns = [c.replace(" ", "_") for c in proteomics_covariates.columns]

    # --- join: keep only samples present in covariates ---
    proteomics_samples = proteomics_covariates.select(["FID", "IID"])
    proteomics_levels = proteomics_levels.join(proteomics_samples, on=["FID", "IID"], how="inner").sort("FID")

    assert proteomics_levels.shape[0] == proteomics_covariates.shape[0], "Sample count mismatch after join"
    assert proteomics_levels["FID"].to_list() == proteomics_covariates["FID"].to_list(), "Sample order mismatch after join"

    # --- write ---
    proteomics_covariates.write_csv(args.output_dir / "regenie_proteomics_covariates.txt", separator="\t")
    proteomics_levels.write_csv(args.output_dir / "regenie_proteomics_levels.txt", separator="\t", null_value="NA")
    proteomics_samples.write_csv(args.output_dir / "regenie_proteomics_samples.txt", separator="\t", include_header=False, null_value="NA")

    print(f"Samples: {proteomics_covariates.shape[0]}, Genes: {len(genes)}, Covariates: {len(covariate_columns)}")


if __name__ == "__main__":
    main()
