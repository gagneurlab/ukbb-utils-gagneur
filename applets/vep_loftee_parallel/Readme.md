Annotates the variant metadata with Ensembl VEP and the LOFTEE plugin. Input variants are split into
chunks that are annotated in parallel subjobs, filtered by gene and transcript biotype, and merged
back into one parquet with one-hot encoded consequences.

The input parquet must have the columns `chrom`, `pos`, `id`, `ref`, `alt` — i.e. the
`variant_metadata_parquet` from `bcf_qc_to_parquet`.

| Input | Type | Default | Description |
|---|---|---|---|
| `variants_parquet` | file | — | Variant metadata to annotate |
| `vep_version` | string | 113 | Ensembl release |
| `chunk_size` | int | 500000 | Variants per subjob |
| `use_loftee` | boolean | true | Run the LOFTEE plugin |
| `use_polyphen` | boolean | true | Run PolyPhen2 predictions |
| `genes_to_keep_file` | file | — | Optional `.txt` of Ensembl gene IDs, one per line |
| `biotypes_filter` | array:string | `["protein_coding"]` | Transcript biotypes to keep; leave empty to keep all |

| Output | Description |
|---|---|
| `vep_parquet` | Variants joined with filtered, one-hot encoded VEP/LOFTEE annotations |