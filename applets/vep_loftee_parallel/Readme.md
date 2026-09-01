<!-- dx-header -->
# VEP and LOFTEE (DNAnexus Platform App)

Annotates variants with Ensembl VEP and the LOFTEE plugin. Input variants are split into chunks that are
annotated in parallel subjobs, filtered by gene and transcript biotype, and merged back into a single
parquet with one-hot encoded consequences.

The input parquet must contain the columns `chrom`, `pos`, `id`, `ref` and `alt` — i.e. the
`variant_metadata_parquet` produced by `bcf_to_gt_parquet`.
<!-- /dx-header -->

Inputs, outputs and example commands: see the [repository README](../../README.md).
