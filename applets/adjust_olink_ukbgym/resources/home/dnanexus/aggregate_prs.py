#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polars",
#   "tqdm",
# ]
# ///

import argparse
from pathlib import Path

import polars as pl
import tqdm


def main():
    parser = argparse.ArgumentParser(description="Aggregate regenie step 1 PRS files into a single parquet.")
    parser.add_argument("--prs-list", required=True, type=Path, help="Path to the .list file (space-separated: gene filename).")
    parser.add_argument("--prs-dir", required=True, type=Path, help="Directory containing the .prs files.")
    parser.add_argument("--output", required=True, type=Path, help="Output parquet path.")
    args = parser.parse_args()

    prslist_df = pl.read_csv(args.prs_list, separator=" ", has_header=False)
    prslist_df.columns = ["gene", "filename"]

    dfs = []
    prev_samples = None
    for row in tqdm.tqdm(prslist_df.iter_rows(named=True), total=prslist_df.height):
        gene = row["gene"]
        filename = row["filename"]
        prs_df = (
            pl.read_csv(args.prs_dir / filename, separator=" ")
            .drop("FID_IID")
            .transpose(header_name="sample", include_header=True, column_names=[gene])
            .drop_nulls()
            .sort("sample")
        )
        curr_samples = prs_df["sample"].to_list()
        if prev_samples is not None:
            assert curr_samples == prev_samples, f"Sample mismatch for gene {gene}"
        prev_samples = curr_samples
        dfs.append(prs_df.drop("sample"))

    prs_df = pl.concat(dfs, how="horizontal").with_columns(
        sample=pl.Series("sample", prev_samples).str.split("_").list.get(0)
    )
    prs_df = prs_df.select(["sample"] + prs_df.columns[:-1])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    prs_df.write_parquet(args.output)
    print(f"Written {prs_df.shape[0]} samples x {prs_df.shape[1]-1} genes to {args.output}")


if __name__ == "__main__":
    main()
