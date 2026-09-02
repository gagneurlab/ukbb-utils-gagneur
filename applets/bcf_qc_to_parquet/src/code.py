#!/usr/bin/env python3

"""
bcf_to_gt_parquet — QC a batch of BCF/VCF files in parallel, then consolidate
into a single gt_long.parquet and variant_metadata.parquet.

Entry points:
  main         — reads input_file_list, fans out process_file subjobs, then
                 waits for all and launches gather
  process_file — QC + normalise one BCF → sparse non-ref GT Parquet
  gather       — concatenate chunk Parquets → gt_long + variant_metadata

TO BUILD:
  dx build bcf_to_gt_parquet/ --destination <PROJECT-ID>:<PATH> -f
"""

import os
import subprocess
import dxpy
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

WORK_DIR = "/home/dnanexus"
LOCAL_BIN = os.path.join(WORK_DIR, ".local", "bin")


# ── Helpers ───────────────────────────────────────────────────────────────────


def run(cmd: str) -> None:
    logger.info("$ %s", cmd)
    env = os.environ.copy()
    env["PATH"] = f"{LOCAL_BIN}:{env['PATH']}"
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, env=env)


def install_tools() -> None:
    """Install bcftools 1.23 and uv (idempotent)."""
    os.makedirs(LOCAL_BIN, exist_ok=True)

    bcftools_bin = os.path.join(LOCAL_BIN, "bcftools")
    if not os.path.exists(bcftools_bin):
        logger.info("Installing bcftools 1.23...")
        run(
            "wget -q https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2 && "
            "tar -xjf bcftools-1.23.tar.bz2 && rm bcftools-1.23.tar.bz2 && "
            "cd bcftools-1.23 && "
            f"./configure --prefix={WORK_DIR}/.local >/dev/null && "
            "make >/dev/null && make install >/dev/null && "
            "cd .. && rm -rf bcftools-1.23"
        )
    else:
        logger.info("bcftools already present, skipping install.")

    uv_bin = os.path.join(LOCAL_BIN, "uv")
    if not os.path.exists(uv_bin):
        logger.info("Installing uv...")
        run(
            "wget -q -O /tmp/uv.tar.gz "
            "https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-unknown-linux-gnu.tar.gz && "
            "tar -xzf /tmp/uv.tar.gz -C /tmp && "
            f"mv /tmp/uv-x86_64-unknown-linux-gnu/uv {LOCAL_BIN}/uv && "
            "rm -rf /tmp/uv.tar.gz /tmp/uv-x86_64-unknown-linux-gnu"
        )
    else:
        logger.info("uv already present, skipping install.")


def strip_bcf_extension(filename: str) -> str:
    for ext in (".bcf.gz", ".vcf.gz", ".bcf"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def _dxid(dxlink):
    """
    Extract (file_id, project_id) from whatever dxpy passes to an entry point.

    dxpy entry points receive file inputs as dxlink dicts, e.g.:
      {"$dnanexus_link": {"id": "file-XXX", "project": "project-YYY"}}
    or the flat form:
      {"$dnanexus_link": "file-XXX"}

    dxpy.download_dxfile / dxpy.DXFile both accept plain string IDs, but
    omit the project context unless it is passed explicitly — which causes
    the /download API to return 404 for files in controlled-access datasets.
    """
    if isinstance(dxlink, dict):
        link = dxlink.get("$dnanexus_link", dxlink)
        if isinstance(link, dict):
            return link["id"], link.get("project")
        return link, None
    return dxlink, None  # already a plain file-XXX string


def dx_name(dxlink) -> str:
    """Return the platform filename for a dxlink."""
    file_id, _ = _dxid(dxlink)
    return dxpy.describe(file_id)["name"]


def dx_download(dxlink, local_path: str) -> None:
    """Download a dxlink to a local path, preserving project context."""
    file_id, project_id = _dxid(dxlink)
    dxpy.download_dxfile(file_id, local_path, project=project_id)


# ── process_file ──────────────────────────────────────────────────────────────


@dxpy.entry_point("process_file")
def process_file(
    input_bcf,
    input_bcf_index,
    af_threshold,
    min_gq,
    min_lad,
    max_missing,
    fasta_ref=None,
    fasta_ref_index=None,
    samples_file=None,
):
    """QC + normalise one BCF/VCF file and write a sparse non-ref GT Parquet."""
    install_tools()

    bcf_fn = dx_name(input_bcf)
    logger.info("Downloading %s ...", bcf_fn)
    dx_download(input_bcf, bcf_fn)

    index_fn = dx_name(input_bcf_index)
    logger.info("Downloading index %s ...", index_fn)
    dx_download(input_bcf_index, index_fn)

    fasta_ref_arg = ""
    if fasta_ref:
        fasta_fn = dx_name(fasta_ref)
        logger.info("Downloading FASTA ref %s ...", fasta_fn)
        dx_download(fasta_ref, fasta_fn)
        if fasta_ref_index:
            fasta_index_fn = dx_name(fasta_ref_index)
            logger.info("Downloading FASTA index %s ...", fasta_index_fn)
            dx_download(fasta_ref_index, fasta_index_fn)
        fasta_ref_arg = f"--fasta-ref {fasta_fn}"

    samples_arg = ""
    if samples_file:
        samples_fn = dx_name(samples_file)
        logger.info("Downloading sample list %s ...", samples_fn)
        dx_download(samples_file, samples_fn)
        with open(samples_fn) as fh:
            n_samples = sum(1 for line in fh if line.strip())
        logger.info("Restricting output to %d samples from %s", n_samples, samples_fn)
        samples_arg = f"--samples-file {samples_fn}"

    output_fn = strip_bcf_extension(bcf_fn) + ".parquet"
    logger.info(
        "Converting %s -> %s  (af=%.4g, min_gq=%s, min_lad=%s, max_missing=%s, samples=%s)",
        bcf_fn,
        output_fn,
        af_threshold,
        min_gq,
        min_lad,
        max_missing,
        samples_fn if samples_file else "all",
    )

    run(
        f"uv run --script {WORK_DIR}/bcf2parquet.py convert "
        f"{bcf_fn} {output_fn} "
        f"--af-threshold {af_threshold} "
        f"--min-gq {min_gq} "
        f"--min-lad {min_lad} "
        f"--max-missing {max_missing} "
        f"{fasta_ref_arg} "
        f"{samples_arg}"
    )

    logger.info("Uploading %s ...", output_fn)
    return {"output_parquet": dxpy.dxlink(dxpy.upload_local_file(output_fn))}


# ── gather ────────────────────────────────────────────────────────────────────


@dxpy.entry_point("gather")
def gather(chunk_parquets):
    """
    Download all chunk Parquets produced by process_file, concatenate into
    gt_long.parquet (id, sample, gt), then derive variant_metadata.parquet
    (id, chrom, pos, ref, alt, sc_cohort, ac_cohort, mac_cohort).
    """
    import polars as pl

    # 1. Download chunk parquets in parallel
    tmp_dir = os.path.join(WORK_DIR, "chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    logger.info("Downloading %d chunk parquets...", len(chunk_parquets))

    def _download(args):
        idx, link = args
        path = os.path.join(tmp_dir, f"chunk_{idx}.parquet")
        try:
            dx_download(link, path)
            return path, True, None
        except Exception as e:
            return path, False, str(e)

    local_paths = []
    with ThreadPoolExecutor(max_workers=16) as ex:
        futures = {
            ex.submit(_download, item): item for item in enumerate(chunk_parquets)
        }
        done = 0
        for fut in as_completed(futures):
            path, ok, err = fut.result()
            done += 1
            if ok:
                local_paths.append(path)
            else:
                logger.warning("Failed to download chunk %s: %s", path, err)
            if done % 200 == 0:
                logger.info("  Downloaded %d/%d chunks...", done, len(chunk_parquets))

    # 2. Validate (skip corrupt files)
    valid_paths = []
    for path in local_paths:
        try:
            pl.scan_parquet(path).head(1).collect()
            valid_paths.append(path)
        except Exception as e:
            logger.warning("Skipping corrupt chunk %s: %s", os.path.basename(path), e)

    logger.info("%d/%d valid chunk parquets.", len(valid_paths), len(chunk_parquets))

    # 3. Consolidate → gt_long.parquet
    #    Each chunk has columns: CHROM, POS, REF, ALT, SAMPLE, GT (from bcf2parquet.py)
    gt_long_path = os.path.join(WORK_DIR, "gt_long.parquet")
    logger.info("Streaming chunks to gt_long.parquet...")

    lfs = [
        pl.scan_parquet(p)
        .with_columns(
            pl.concat_str(
                [
                    pl.col("CHROM"),
                    pl.lit(":"),
                    pl.col("POS").cast(pl.Utf8),
                    pl.lit(":"),
                    pl.col("REF"),
                    pl.lit(":"),
                    pl.col("ALT"),
                ]
            ).alias("id")
        )
        .select(["id", "SAMPLE", "GT"])
        .rename({"SAMPLE": "sample", "GT": "gt"})
        for p in valid_paths
    ]
    pl.concat(lfs, how="vertical_relaxed").sink_parquet(
        gt_long_path, engine="streaming"
    )
    logger.info("gt_long.parquet written.")

    # 4. Compute per-variant counts → variant_metadata.parquet
    tmp_meta_path = os.path.join(WORK_DIR, "tmp_variant_metadata.parquet")
    meta_path = os.path.join(WORK_DIR, "variant_metadata.parquet")

    logger.info("Computing per-variant sample count and allele count...")
    (
        pl.scan_parquet(gt_long_path)
        .select(["id", "gt"])
        .group_by("id")
        .agg(
            sc_cohort=pl.len(),
            ac_cohort=pl.sum("gt"),
        )
        .sink_parquet(tmp_meta_path, engine="streaming")
    )

    max_ac = pl.scan_parquet(tmp_meta_path).select("ac_cohort").max().collect().item()
    logger.info("Max AC in cohort: %d. Computing MAC and splitting id...", max_ac)

    (
        pl.scan_parquet(tmp_meta_path)
        .with_columns(
            chrom=pl.col("id").str.split(":").list.get(0),
            pos=pl.col("id").str.split(":").list.get(1).cast(pl.Int64),
            ref=pl.col("id").str.split(":").list.get(2),
            alt=pl.col("id").str.split(":").list.get(3),
            mac_cohort=(
                pl.when(pl.col("ac_cohort") > int(max_ac / 2))
                .then(max_ac + 1 - pl.col("ac_cohort"))
                .otherwise(pl.col("ac_cohort"))
            ),
        )
        .sink_parquet(meta_path, engine="streaming")
    )
    logger.info("variant_metadata.parquet written.")

    # 5. Upload outputs
    logger.info("Uploading outputs...")
    return {
        "gt_long_parquet": dxpy.dxlink(dxpy.upload_local_file(gt_long_path)),
        "variant_metadata_parquet": dxpy.dxlink(dxpy.upload_local_file(meta_path)),
    }


# ── main ──────────────────────────────────────────────────────────────────────


@dxpy.entry_point("main")
def main(
    input_file_list,
    af_threshold=0.001,
    min_gq=10,
    min_lad=8,
    max_missing=0.1,
    fasta_ref=None,
    fasta_ref_index=None,
    samples_file=None,
):
    """
    Read a CSV or Parquet file list and dispatch one process_file subjob per row,
    then a single gather subjob that waits on all of them.

    Required columns in input_file_list: vcf_file_id, vcf_index_id
    Values must be DNAnexus file IDs (e.g. file-GZZ3630J55kk...).

    samples_file, if given, restricts the output cohort — see process_file. Note
    that the sc_/ac_/mac_cohort columns of variant_metadata.parquet are then
    computed over that cohort, since gather derives them from gt_long.
    """
    import polars as pl

    list_path = os.path.join(WORK_DIR, "input_file_list")
    list_name = dx_name(input_file_list)
    dx_download(input_file_list, list_path)

    df = (
        pl.read_parquet(list_path)
        if list_name.endswith(".parquet")
        else pl.read_csv(list_path)
    )

    missing_cols = {"vcf_file_id", "vcf_index_id"} - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"input_file_list is missing required columns: {missing_cols}. Found: {df.columns}"
        )

    def _link(file_ref: str):
        """Parse 'project-XXX:file-YYY' or bare 'file-YYY' into a valid dxlink.

        Bare file IDs are assumed to live in the current job's project
        (dxpy.WORKSPACE_ID), which ensures the project context is always
        included in the dxlink and avoids 404s on controlled-access data.
        """
        if ":" in file_ref:
            project_id, file_id = file_ref.split(":", 1)
        else:
            project_id = dxpy.PROJECT_CONTEXT_ID
            file_id = file_ref
        return dxpy.dxlink(file_id, project_id=project_id)

    logger.info("Dispatching %d process_file subjobs...", len(df))

    chunk_outputs = []
    for row in df.iter_rows(named=True):
        fn_input = {
            "input_bcf": _link(row["vcf_file_id"]),
            "input_bcf_index": _link(row["vcf_index_id"]),
            "af_threshold": af_threshold,
            "min_gq": min_gq,
            "min_lad": min_lad,
            "max_missing": max_missing,
        }
        if fasta_ref is not None:
            fn_input["fasta_ref"] = fasta_ref
        if fasta_ref_index is not None:
            fn_input["fasta_ref_index"] = fasta_ref_index
        if samples_file is not None:
            fn_input["samples_file"] = samples_file

        subjob = dxpy.new_dxjob(
            fn_name="process_file",
            fn_input=fn_input,
            instance_type={"*": "mem2_ssd1_v2_x2", "2": "mem3_ssd1_v2_x2"},
        )
        chunk_outputs.append(subjob.get_output_ref("output_parquet"))

    logger.info(
        "Launching gather subjob (depends on all %d process_file jobs)...",
        len(chunk_outputs),
    )
    gather_job = dxpy.new_dxjob(
        fn_name="gather",
        fn_input={"chunk_parquets": chunk_outputs},
        instance_type="mem2_ssd2_v2_x32",
    )

    return {
        "gt_long_parquet": gather_job.get_output_ref("gt_long_parquet"),
        "variant_metadata_parquet": gather_job.get_output_ref(
            "variant_metadata_parquet"
        ),
    }


dxpy.run()
