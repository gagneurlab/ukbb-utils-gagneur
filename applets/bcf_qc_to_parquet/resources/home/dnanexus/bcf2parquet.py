#!/usr/bin/env -S uv run --script

# /// script
# requires-python = '>=3.11'
# dependencies = [
#   "polars>1.3",
#   "click>=8.0",
# ]
# ///

import logging
import polars as pl
import subprocess
import click
import tempfile
from polars.exceptions import NoDataError  # still used to detect empty TSV

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Schema for fast parsing
SCHEMA = {
    "CHROM": pl.Utf8,
    "POS": pl.Int64,
    "REF": pl.Utf8,
    "ALT": pl.Utf8,
    "SAMPLE": pl.Utf8,
    "GT": pl.Utf8,  # Read as string first to handle the / and |
}


def bcf_to_tsv(
    bcf_file: str,
    af_threshold: float,
    output_file,
    min_gq: int | None = None,
    min_lad: int | None = None,
    max_missing: float | None = None,
    fasta_ref: str | None = None,
    samples_file: str | None = None,
) -> None:
    """
    Process BCF/VCF file through bcftools pipeline and return TSV data.

    Pipeline: QC (annotate + setGT + missingness) -> normalize -> MAF filter
              -> cohort subset -> sparse GT output
    """
    steps = []

    # 1. Keep only relevant FORMAT fields (GT, GQ, LAD)
    steps.append(
        f"bcftools annotate -x ^FORMAT/GT,^FORMAT/GQ,^FORMAT/LAD -Ou {bcf_file}"
    )

    # 2. Set low-quality genotypes to missing
    setgt_conditions = []
    if min_gq is not None:
        setgt_conditions.append(f"FMT/GQ<={min_gq}")
    if min_lad is not None:
        setgt_conditions.append(f"smpl_sum(FMT/LAD)<{min_lad}")

    if setgt_conditions:
        filter_expr = " | ".join(setgt_conditions)
        steps.append(f'bcftools +setGT -Ou -- -t q -i "{filter_expr}" -n .')

    # 3. Filter variants with too much missingness
    if max_missing is not None:
        steps.append(f'bcftools filter -e "F_MISSING > {max_missing}" -Ou')

    # 4. Normalize and split multiallelics
    if fasta_ref is not None:
        steps.append(f"bcftools norm -m - --fasta-ref {fasta_ref} -Ou")
    else:
        steps.append("bcftools norm -m - -Ou")

    # 5. MAF filter
    steps.append("bcftools +fill-tags -Ou -- -t AF,MAF")
    steps.append(f'bcftools view -i "MAF < {af_threshold}" -Ou')

    # 6. Restrict to the analysis cohort.
    #
    #    Deliberately placed AFTER the MAF and missingness filters: those tags are
    #    computed over whatever samples are in the stream, so subsetting earlier
    #    would define the variant set on the subset rather than on the full cohort.
    #    Keeping it here reproduces the variant set of the published run while
    #    still avoiding TSV rows for samples that get discarded anyway.
    #
    #    No --force-samples: an ID in the list that is absent from the BCF should
    #    fail the job loudly rather than silently yield a smaller cohort.
    if samples_file is not None:
        steps.append(f"bcftools view --samples-file {samples_file} -Ou")

    # 7. Extract sparse non-ref genotypes
    steps.append(
        "bcftools query "
        '--include \'GT!="RR" & GT!="mis" & GT!="R"\' '
        "--format '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]'"
    )

    cmd = "set -o pipefail; " + " | ".join(steps)
    logger.info("bcftools pipeline:\n  %s", "\n  | ".join(steps))
    subprocess.run(
        cmd, shell=True, executable="/bin/bash", stdout=output_file, check=True
    )
    logger.info("bcftools pipeline complete")


def tsv_to_parquet(tsv_file, output_file: str):
    """
    Convert TSV data to Parquet format.

    Args:
        tsv_data: TSV data as bytes
        output_file: Output Parquet file path
    """
    logger.info("Converting TSV to Parquet...")
    try:
        df = pl.scan_csv(
            tsv_file,
            separator="\t",
            has_header=False,
            new_columns=["CHROM", "POS", "REF", "ALT", "SAMPLE", "GT"],
            infer_schema=False,
            schema=SCHEMA,
        ).with_columns(
            # Replace VCF GT strings with integers and cast to Int8
            pl.col("GT")
            # 0/1, 1/0, 0|1, 1|0 -> "1"
            .str.replace(r"(0[|/]1)|(1[|/]0)", "1")
            # 1/1 or 1|1 -> "2"
            .str.replace(r"1[|/]1", "2")
            .cast(pl.Int8, strict=False)
        )
    except NoDataError:
        raise RuntimeError(
            f"TSV output is empty — bcftools pipeline may have been OOM-killed "
            f"or all variants were filtered out. Check {tsv_file} and job logs."
        )

    # Save to Parquet
    df.sink_parquet(output_file, compression="zstd")


@click.group()
def cli():
    """BCF to Parquet converter tool."""
    pass


@cli.command()
def install():
    """
    Install dependencies (no-op when using pre-installed packages).
    """
    click.echo("Dependencies installed")


@cli.command()
@click.argument("bcf_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option(
    "--af-threshold",
    default=0.001,
    type=float,
    help="Allele frequency threshold for filtering (default: 0.001)",
)
@click.option(
    "--min-gq",
    default=10,
    type=int,
    help="Min GQ; at or below is set to missing (default: 10)",
)
@click.option(
    "--min-lad",
    default=8,
    type=int,
    help="Min sum of LAD allele depths; below is set to missing (default: 8)",
)
@click.option(
    "--max-missing",
    default=0.1,
    type=float,
    help="Max fraction of missing genotypes per variant (default: 0.1)",
)
@click.option(
    "--fasta-ref",
    default=None,
    type=click.Path(exists=True),
    help="Reference FASTA for bcftools norm",
)
@click.option(
    "--samples-file",
    default=None,
    type=click.Path(exists=True),
    help=(
        "Restrict output to these samples (one ID per line, no header). "
        "Applied after the MAF/missingness filters, so the variant set stays "
        "defined on the full cohort."
    ),
)
def convert(
    bcf_file: str,
    output_file: str,
    af_threshold: float,
    min_gq: int,
    min_lad: int,
    max_missing: float,
    fasta_ref: str | None,
    samples_file: str | None,
):
    """
    Convert BCF/VCF file to Parquet format.

    Applies QC (FORMAT field pruning, GQ/LAD filtering, missingness), normalization,
    MAF filtering, an optional cohort subset, then outputs sparse non-ref genotypes
    as Parquet.

    \b
    Example:
        python bcf2parquet.py convert input.vcf.gz output.parquet \\
            --min-gq 10 --min-lad 8 --max-missing 0.1 \\
            --fasta-ref ref.fa --af-threshold 0.001 \\
            --samples-file eur_samples.txt
    """
    logger.info("Processing %s", bcf_file)
    logger.info(
        "Parameters: af_threshold=%.4g, min_gq=%s, min_lad=%s, max_missing=%s, "
        "fasta_ref=%s, samples_file=%s",
        af_threshold,
        min_gq,
        min_lad,
        max_missing,
        fasta_ref or "none",
        samples_file or "none",
    )

    with tempfile.NamedTemporaryFile(mode="w+b") as tmp:
        bcf_to_tsv(
            bcf_file,
            af_threshold,
            tmp,
            min_gq=min_gq,
            min_lad=min_lad,
            max_missing=max_missing,
            fasta_ref=fasta_ref,
            samples_file=samples_file,
        )
        result = subprocess.run(["du", "-sh", tmp.name], capture_output=True, text=True)
        logger.info("TSV temp file size: %s", result.stdout.strip())
        tsv_to_parquet(tmp.name, output_file)

    logger.info("Done -> %s", output_file)


if __name__ == "__main__":
    cli()
