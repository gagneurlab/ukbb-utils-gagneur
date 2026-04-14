# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polars",
#   "mygene",
# ]
# ///
"""Standalone script: preprocess raw Olink proteomics data to Ensembl-column parquet.

Usage:
    python preprocess_olink.py \
        --olink-path       <path/to/olink.csv> \
        --helper-assay-path <path/to/olink_helper_assay.tsv> \
        --output-dir       <path/to/output_dir> \
        [--samples-path    <path/to/samples.csv>]

Outputs (written to --output-dir):
    olink_preprocessed.parquet   – wide parquet, columns: sample + Ensembl IDs
    olink_genes.txt              – one Ensembl ID per line
    olink_samples.txt            – one sample ID per line
"""

import argparse
import logging
from collections import Counter
from pathlib import Path

import polars as pl
from mygene import MyGeneInfo

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_valid_chr(chr: str) -> bool:
    valid_chrs = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    return chr in valid_chrs


def _extract_gene_from_result(result: dict) -> str | list | None:
    ensembl = result.get("ensembl")
    if isinstance(ensembl, dict):
        return ensembl.get("gene")
    elif isinstance(ensembl, list):
        genes = [g.get("gene") for g in ensembl]
        genomic_pos = result.get("genomic_pos", [])
        valid_genes = [g for i, g in enumerate(genes) if _is_valid_chr(genomic_pos[i].get("chr"))]
        if len(valid_genes) == 1:
            return valid_genes[0]
        elif len(valid_genes) == 0:
            return None
        else:
            return valid_genes
    return None


def load_assay_uniprot_mapping(assay_file_path: Path) -> dict[str, str]:
    df = pl.read_csv(assay_file_path, separator="\t", has_header=True)
    mapping = {row[0].lower().replace("-", "_"): row[1] for row in df.rows()}

    uniprot_counts = Counter(mapping.values())
    multi = {uid: c for uid, c in uniprot_counts.items() if c > 1}
    if multi:
        logger.info("UniProt IDs with more than one assay: %d", len(multi))
    return mapping


def map_uniprot_to_ensembl(uniprot_ids: set[str]) -> tuple[dict[str, str], set[str], set[str]]:
    mg = MyGeneInfo()
    logger.info("Querying %d UniProt IDs via MyGeneInfo...", len(uniprot_ids))
    results = mg.querymany(
        uniprot_ids,
        scopes="uniprot",
        fields=["ensembl.gene", "genomic_pos.chr"],
        species="human",
    )

    mapping: dict[str, str] = {}
    dup_ids: set[str] = set()
    missing_ids: set[str] = set()

    for result in results:
        uid = result["query"].split("-")[0]
        gene = _extract_gene_from_result(result)

        if uid in mapping:
            dup_ids.add(uid)
            mapping.pop(uid)
            continue

        if gene is None:
            if uid not in mapping:
                missing_ids.add(uid)
        elif isinstance(gene, str) and uid not in dup_ids:
            mapping[uid] = gene
            missing_ids.discard(uid)
        elif isinstance(gene, list):
            dup_ids.add(uid)

    unmapped = len(uniprot_ids) - len(mapping)
    logger.info("Mapped: %d / %d  (missing=%d, duplicate=%d)",
                len(mapping), len(uniprot_ids), len(missing_ids), len(dup_ids))
    assert unmapped == len(missing_ids) + len(dup_ids)
    return mapping, dup_ids, missing_ids


def create_assay_ensembl_mapping(
    assay_uniprot: dict[str, str],
    uniprot_ensembl: dict[str, str],
) -> dict[str, str]:
    return {
        assay: uniprot_ensembl[uid]
        for assay, uid in assay_uniprot.items()
        if uid in uniprot_ensembl
    }


def load_and_transform_olink(
    olink_path: Path,
    assay_ensembl: dict[str, str],
    assay_uniprot: dict[str, str],
) -> tuple[pl.LazyFrame, list[str], list[str]]:
    logger.info("Loading Olink data from %s", olink_path)
    df = pl.scan_csv(olink_path).rename({"eid": "sample"}).with_columns(pl.col("sample").cast(pl.String))

    assays_in_data = [c for c in df.columns if c != "sample"]
    missing = [a for a in assay_uniprot if a not in assays_in_data]
    invalid = [a for a in assays_in_data if a not in assay_uniprot]
    logger.info("Assays missing from data: %d  |  invalid in data: %d", len(missing), len(invalid))

    df = df.select("sample", *assay_ensembl.keys()).rename(assay_ensembl)
    gene_cols = [c for c in df.columns if c != "sample"]
    samples = df.select("sample").collect()["sample"].to_list()
    logger.info("Transformed: %d samples, %d genes", len(samples), len(gene_cols))
    return df, gene_cols, samples


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--olink-path", required=True, type=Path, help="Path to raw Olink CSV file")
    parser.add_argument("--helper-assay-path", required=True, type=Path, help="Path to Olink helper assay TSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to write outputs")
    parser.add_argument("--samples-path", type=Path, default=None, help="Optional CSV of sample IDs to subset to")
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    # Step 1 – assay → UniProt
    assay_uniprot = load_assay_uniprot_mapping(args.helper_assay_path)

    # Step 2 – UniProt → Ensembl
    uniprot_ensembl, _, _ = map_uniprot_to_ensembl(set(assay_uniprot.values()))

    # Step 3 – assay → Ensembl
    assay_ensembl = create_assay_ensembl_mapping(assay_uniprot, uniprot_ensembl)

    # Step 4 – load & transform
    olink_df, ensembl_ids, olink_samples = load_and_transform_olink(
        args.olink_path, assay_ensembl, assay_uniprot
    )

    # Optional sample subsetting
    if args.samples_path:
        logger.info("Subsetting to samples in %s", args.samples_path)
        samples = list(
            pl.read_csv(args.samples_path, has_header=True)
            .select(pl.nth(0).cast(pl.String))
            .to_series()
        )
        missing_samples = set(samples) - set(olink_samples)
        if missing_samples:
            logger.warning("%d samples not found in Olink data and will be skipped: %s", len(missing_samples), missing_samples)
            samples = [s for s in samples if s not in missing_samples]
        sample_df = pl.DataFrame({"sample": samples})
        olink_df = sample_df.join(olink_df.collect(), on="sample", how="left")
        olink_samples = samples

    # Write outputs
    olink_out = out / "olink_preprocessed.parquet"
    genes_out = out / "olink_genes.txt"
    samples_out = out / "olink_samples.txt"

    if isinstance(olink_df, pl.LazyFrame):
        olink_df = olink_df.collect()
    olink_df.write_parquet(olink_out)
    logger.info("Saved parquet: %s", olink_out)

    genes_out.write_text("\n".join(ensembl_ids) + "\n")
    logger.info("Saved gene list: %s", genes_out)

    samples_out.write_text("\n".join(str(s) for s in olink_samples) + "\n")
    logger.info("Saved sample list: %s", samples_out)


if __name__ == "__main__":
    main()
