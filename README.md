# DeepRVAT-WGS

DNAnexus applets for running [DeepRVAT](https://github.com/PMBio/deeprvat) on UK Biobank whole-genome
sequencing data on the [UKB Research Analysis Platform (RAP)](https://ukbiobank.dnanexus.com/).

Each applet is one step of the pipeline and runs on the RAP as a standalone job.

```
BCF/VCF files ──▶ vcf_qc_to_parquet ──▶ vep_loftee_parallel ──▶ annotated variants
                                                                        +
UKB dataset  ──▶ extract_phenotypes_and_covariates ─────────────▶ association testing

Olink (optional): prepare_regenie_olink_inputs ──▶ regenie ──▶ adjust_olink_ukbgym
```

## Setup

```bash
pip install dxpy pre-commit
dx login          # authenticate against the RAP
dx select         # pick a project, e.g. ukb-gagneur
pre-commit install
```

Pre-commit strips notebook outputs and applies formatting on every commit — please install it.

## Building and running an applet

Build once (from the repo root), then run as often as needed:

```bash
# build / update the applet on the platform
dx build -f applets/<applet_dir>/ --destination <project-id>:/applets/

# run it
dx run <applet_name> -i<input>=<value> --destination /path/for/outputs/
```

`<applet_dir>` is the folder name below; `<applet_name>` is the `name` field in its `dxapp.json`
(the same as the folder name, except `vcf_qc_to_parquet`, which builds as `bcf_to_gt_parquet`).

Applets can also be run from the RAP web UI: select the applet, fill in the inputs, click **Run**.

---

## `vcf_qc_to_parquet` — QC genotypes and convert to parquet

Runs genotype- and variant-level QC on a batch of BCF/VCF files (one parallel subjob per file) and
consolidates the results into two parquet tables. Filters on genotype quality, allele depth, per-variant
missingness and MAF, and normalises against a reference FASTA if one is given.

Built and run as **`bcf_to_gt_parquet`**.

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

```bash
dx run bcf_to_gt_parquet \
  -iinput_file_list="/path/to/bcf_file_list.csv" \
  -iaf_threshold=0.001 \
  --destination /processed_data/genotypes/
```

---

## `vep_loftee_parallel` — annotate variants with VEP + LOFTEE

Annotates the variant metadata with Ensembl VEP and the LOFTEE plugin. Input variants are split into
chunks that are annotated in parallel subjobs, filtered by gene and transcript biotype, and merged
back into one parquet with one-hot encoded consequences.

The input parquet must have the columns `chrom`, `pos`, `id`, `ref`, `alt` — i.e. the
`variant_metadata_parquet` from `vcf_qc_to_parquet`.

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

```bash
dx run vep_loftee_parallel \
  -ivariants_parquet="/processed_data/genotypes/variant_metadata.parquet" \
  -ivep_version="113" \
  -ichunk_size=500000 \
  -iuse_loftee=true \
  --destination /processed_data/annotations/
```

---

## `extract_phenotypes_and_covariates` — phenotypes, covariates and PRS

Extracts phenotype and covariate fields straight from the UKB Apollo dataset (batched to stay under the
API limits), optionally subsets to European ancestry and unrelated individuals, and applies statin
correction to cholesterol/LDL. If PLINK files are supplied it also runs REGENIE and returns PRS and
PRS-corrected phenotypes.

Field lists for the UKBGym analysis live in the applet directory —
`ukbbgym_trait_fieldIDs.txt` and `ukbbgym_covariate_fieldIDs.txt`. These are repository files, not part
of the built applet, so upload them to the project first and pass the resulting file IDs.

| Input | Type | Default | Description |
|---|---|---|---|
| `dataset_id` | record | — | UKB dataset record ID (`record-XXXX`) |
| `pheno_list_file` | file | — | Text file listing phenotype columns |
| `additional_covars_file` | file | — | Optional text file listing extra covariate columns |
| `bed_file` / `bim_file` / `fam_file` | file | — | Optional PLINK files; without them no PRS or corrected phenotypes are produced |
| `subset_eur` | boolean | true | Keep European ancestry only (`p30079 = 1`) |
| `subset_unrelated` | boolean | true | Drop related samples via the kinship matrix |
| `statin_correction` | boolean | true | Apply statin correction to lipid phenotypes |

| Output | Description |
|---|---|
| `phenotypes`, `covariates` | Extracted and filtered tables (parquet) |
| `PRS`, `corrected_phenotypes` | REGENIE PRS and PRS-corrected phenotypes (only if PLINK files were given) |

```bash
# upload the field lists once; reuse the file IDs on later runs
pheno_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_trait_fieldIDs.txt --brief)
covar_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_covariate_fieldIDs.txt --brief)

dx run extract_phenotypes_and_covariates \
  -idataset_id="record-XXXXXXXX" \
  -ipheno_list_file="$pheno_file" \
  -iadditional_covars_file="$covar_file" \
  -ibed_file="/plink/ukb.bed" -ibim_file="/plink/ukb.bim" -ifam_file="/plink/ukb.fam" \
  --destination /processed_data/phenotypes/
```

---

## `prepare_regenie_olink_inputs` — format Olink data for REGENIE

Maps raw UKB covariate fields to readable labels, one-hot encodes categoricals, adds derived features
(age², age×sex, …), joins them to the Olink protein levels and writes the three tab-separated files
REGENIE expects.

Export the inputs from the dataset first with `table-exporter`:

```bash
fields_file=$(dx upload olink_covariate_fields.txt --brief)

dx run table-exporter \
  -idataset_or_cohort_or_dashboard=record-REDACTED \
  -ientity="participant" \
  -ifield_names_file_txt="$fields_file" \
  -ioutput="olink_covariates_export" \
  -icoding_option="RAW" -iheader_style="FIELD-NAME" \
  -
  --destination="/processed_data/olink/raw/"
```

| Input | Type | Description |
|---|---|---|
| `raw_covariates` | file | Raw covariates CSV from `table-exporter` (samples × UKB fields) |
| `covariate_mapping` | file | `olink_covariate_field_mapping.parquet` — maps UKB field IDs to names/titles |
| `proteomics_levels` | file | Olink protein levels CSV from `table-exporter` |

| Output | Description |
|---|---|
| `regenie_covariates` | FID, IID, covariates |
| `regenie_levels` | FID, IID, proteins |
| `regenie_samples` | FID, IID |

```bash
dx run prepare_regenie_olink_inputs \
  -iraw_covariates="/processed_data/olink/raw/olink_covariates_export.csv" \
  -icovariate_mapping="/processed_data/olink/olink_covariate_field_mapping.parquet" \
  -iproteomics_levels="/processed_data/olink/raw/olink_export.csv" \
  --destination /processed_data/olink/regenie/
```

Run REGENIE step 1 on these files to produce the `.prs` files needed by the next applet.

---

## `adjust_olink_ukbgym` — PRS-adjust Olink proteomics

Normalises the Olink protein levels, regresses out the REGENIE polygenic risk scores per protein (OLS),
and runs [Protrider](https://github.com/gtsitsiridis/protadjust) — an autoencoder — on the residuals to
flag outliers. The applet is a thin wrapper around the `ghcr.io/gtsitsiridis/protadjust` Docker image,
pulled fresh on every run, so a new image release needs no rebuild.

| Input | Type | Description |
|---|---|---|
| `input_csv` | array:file | Olink CSV(s) exported with `table-exporter` (entity `olink_instance_0`) |
| `regenie_step_1_prs` | array:file | The `.prs` files from REGENIE step 1 |
| `regenie_step_1_prs_list` | file | The `*_prs.list` file mapping proteins to `.prs` filenames |
| `sample_list` | file | Line-separated sample IDs to include |

| Output | Description |
|---|---|
| `output_parquet` | `adjusted_proteomics.parquet` — PRS-adjusted protein abundances |

```bash
dx run adjust_olink_ukbgym \
  -iinput_csv="/processed_data/olink/raw/olink_export.csv" \
  -iregenie_step_1_prs_list="/regenie/step1_prs.list" \
  -isample_list="/processed_data/olink/samples.txt" \
  --destination /processed_data/olink/adjusted/
```

Protrider is memory-hungry; if the default `mem1_ssd1_v2_x16` is not enough, override it at launch with
`--instance-type mem2_ssd1_v2_x32`.

---

## Tips

- Everything on a RAP workstation or Jupyter session is **deleted when the job ends** — upload results
  back to the project before the session expires.
- `dx watch <job-id>` follows a running job's log; `dx describe <job-id>` shows its inputs and cost.

## Citation

If you use this software, please cite it using the metadata in [CITATION.cff](CITATION.cff).
