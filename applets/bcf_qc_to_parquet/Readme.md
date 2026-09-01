# BCF to GT Parquet (DNAnexus Platform App)

Runs genotype- and variant-level QC on a batch of BCF/VCF files (one parallel subjob per file) and
consolidates the results into a long parquet file. Filters on genotype quality, allele depth, per-variant
missingness and MAF, and normalises against a reference FASTA if one is given.


| Input | Type | Default | Description |
|---|---|---|---|
| `input_file_list` | file | — | CSV/parquet with columns `vcf_file_id`, `vcf_index_id`; one row per BCF + index pair |
| `af_threshold` | float | 0.001 | MAF threshold |
| `min_gq` | int | 10 | Genotypes at or below this GQ are set to missing |
| `min_lad` | int | 8 | Genotypes with LAD sum below this are set to missing |
| `max_missing` | float | 0.1 | Max fraction of missing genotypes per variant |
| `fasta_ref` / `fasta_ref_index` | file | — | Optional reference FASTA + `.fai` for `bcftools norm` |

| Output | Description |
|---|---|
| `gt_long_parquet` | Sparse non-ref genotypes: `id`, `sample`, `gt` |
| `variant_metadata_parquet` | Per-variant stats: `id`, `chrom`, `pos`, `ref`, `alt`, `sc_cohort`, `ac_cohort`, `mac_cohort` |