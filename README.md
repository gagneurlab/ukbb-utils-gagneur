# DeepRVAT-WGS

DNAnexus applets for working with the UKB Biobank data.

Each applet is one step of the pipeline and runs on the RAP as a standalone job.

```
BCF/VCF files ──▶ vcf_qc_to_parquet ──▶ vep_loftee_parallel ──▶ annotated variants
                                                                        +
UKB dataset  ──▶ extract_phenotypes_and_covariates ─────────────▶ association testing

Olink: prepare_regenie_olink_inputs ──▶ regenie ──▶ adjust_olink_ukbgym
```

## Setup

```bash
dx login          # authenticate against the RAP
dx select         # pick a project, e.g. ukb-gagneur
```

## Building and running an applet

Build once (from the repo root), then run as often as needed:

```bash
# build / update the applet on the platform
dx build -f applets/<applet_dir>/ --destination <project-id>:/applets/

# run it
dx run <applet_name> -i<input>=<value> --destination /path/for/outputs/
```

`<applet_dir>` is the folder name below; `<applet_name>` is the `name` field in its `dxapp.json`

Applets can also be run from the RAP web UI.

---

## `bcf_qc_to_parquet` — QC BCF WGS genotypes and convert to parquet
<!-- TODO add the file list for ukbgym -->
Runs genotype- and variant-level QC on a batch of BCF/VCF files and
consolidates the results into a long genotype parquet.

```bash
base=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
wget $base/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget $base/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

dx mkdir -p /reference/
dx upload GRCh38_full_analysis_set_plus_decoy_hla.fa     --destination /reference/
dx upload GRCh38_full_analysis_set_plus_decoy_hla.fa.fai --destination /reference/
```

```bash
dx run bcf_to_gt_parquet \
  -iinput_file_list="/path/to/bcf_file_list.csv" \
  -iaf_threshold=0.001 \
  -ifasta_ref="/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  -ifasta_ref_index="/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" \
  --destination /processed_data/genotypes/
```

---

## `vep_loftee_parallel` — annotate variants with VEP + LOFTEE

Annotates the variant metadata with Ensembl VEP and the LOFTEE plugin. Input variants are split into
chunks that are annotated in parallel subjobs.


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

Extracts phenotype and covariate fields from the UKB dataset, subsets to European ancestry and unrelated individuals, and applies statin correction to cholesterol/LDL. 
This also runs REGENIE Step1 to obtain polygenic risk scores which are then regressed out from the phenotypes, along with additional covariates listed in `ukbbgym_covariate_fieldIDs.txt`.

```bash
# upload the field lists once; reuse the file IDs on later runs
pheno_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_trait_fieldIDs.txt --brief)
covar_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_covariate_fieldIDs.txt --brief)

dx run extract_phenotypes_and_covariates \
  -idataset_id="record-XXXXXXXX" \
  -ipheno_list_file="$pheno_file" \
  -iadditional_covars_file="$covar_file" \
  -ibed_file="/plink/ukb.bed" -ibim_file="/plink/ukb.bim" -ifam_file="/plink/ukb.fam" \
  -isubset_eur="True" \
  -isubset_unrelated="True" \
  -istatin_correction="True" \
  --destination /processed_data/phenotypes/
```

---

## `prepare_regenie_olink_inputs` — format Olink data for REGENIE
<!-- TODO add the file list for olink -->

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
