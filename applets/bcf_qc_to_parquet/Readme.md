<!-- dx-header -->
# BCF to GT Parquet (DNAnexus Platform App)

Runs genotype- and variant-level QC on a batch of BCF/VCF files — one parallel subjob per file — and
consolidates the results into two parquet tables: sparse non-ref genotypes and per-variant metadata.

Filters on genotype quality (`min_gq`), allele depth (`min_lad`), per-variant missingness
(`max_missing`) and minor allele frequency (`af_threshold`). If a reference FASTA is supplied, variants
are normalised with `bcftools norm` first.

Builds and runs as **`bcf_to_gt_parquet`**.
<!-- /dx-header -->

Inputs, outputs and example commands: see the [repository README](../../README.md).
