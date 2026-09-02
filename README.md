# UKBGym processing applets

DNAnexus applets for working with UKB data, particularly for running the UKBB analyses of [UKBGym](https://github.com/gagneurlab/ukbbgym/tree/main/ukbb). 
Each applet is one step of the pipeline and runs on the RAP as a standalone job.
All applets listed below need to be run in order to prepare required UKBGym data. 


```
Phenotypes and Covariates ──▶ extract_phenotypes_and_covariates ──┬──▶ corrected phenotypes
                                                                  │
                                                     eur_samples.txt
                                                                  │
                                                                  ▼
Genotype BCF/VCF files ─────────────────────▶ bcf_qc_to_parquet ──▶ vep_loftee_parallel ──▶ annotated variants

Olink: prepare_regenie_olink_inputs ──▶ regenie ──▶ adjust_olink_ukbgym
```

## Setup

```bash
pip install dxpy
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

## Input file lists

Every list-type input the applets need — file manifests, UKB field IDs, gene and sample lists — lives
in [`ukbgym_file_lists/`](ukbgym_file_lists/).Upload the whole directory to the project once and point the applets at it from there:

```bash
dx upload -r ukbgym_file_lists/ --destination /
```

That creates `/ukbgym_file_lists/` in the project, which is the path used in every example below.
| File | Used by | As input |
|---|---|---|
| `ukbgym_gene_vcf_files.parquet` | `bcf_qc_to_parquet` | `input_file_list` |
| `ukbbgym_gene_list.txt` | `vep_loftee_parallel` | `genes_to_keep_file` |
| `ukbbgym_trait_fieldIDs.txt` | `extract_phenotypes_and_covariates` | `pheno_list_file` |
| `ukbbgym_covariate_fieldIDs.txt` | `extract_phenotypes_and_covariates` | `additional_covars_file` |
| `olink_covariate_fields.txt` | `table-exporter` (Olink covariates) | `field_names_file_txt` |
| `olink_fields.txt` | `table-exporter` (Olink levels) | `field_names_file_txt` |
| `olink_covariate_field_mapping.parquet` | `prepare_regenie_olink_inputs` | `covariate_mapping` |
| `olink_samples.txt` | `adjust_olink_ukbgym` | `sample_list` |


---

## Extract phenotypes and covariates 

### `extract_phenotypes_and_covariates` — phenotypes, covariates and PRS

Extracts phenotype and covariate fields from the UKB dataset, subsets to European ancestry and unrelated individuals, and applies statin correction to cholesterol/LDL. 
This also runs REGENIE Step1 to obtain polygenic risk scores which are then regressed out from the phenotypes, along with additional covariates listed in `ukbbgym_covariate_fieldIDs.txt`.
Needs no input from any other step, so run it first — it emits the sample lists the later steps
take as input:

| Output | Cohort | Feed it to |
|---|---|---|
| `eur_samples.txt` | after the ancestry filter, **before** relatedness pruning | `bcf_qc_to_parquet` → `samples_file` |
| `eur_unrelated_samples.txt` | after ancestry filtering **and** relatedness pruning | `avg_pheno_per_variant_traits` → `samples_txt` |


```bash
dx run extract_phenotypes_and_covariates \
  -idataset_id="record-XXXXXXXX" \
  -ipheno_list_file="/ukbgym_file_lists/ukbbgym_trait_fieldIDs.txt" \
  -iadditional_covars_file="/ukbgym_file_lists/ukbbgym_covariate_fieldIDs.txt" \
  -ibed_file="/plink/ukb.bed" -ibim_file="/plink/ukb.bim" -ifam_file="/plink/ukb.fam" \
  -isubset_eur="True" \
  -isubset_unrelated="True" \
  -istatin_correction="True" \
  --destination /processed_data/phenotypes/
```

---

## Prepare and annotate genotypes

### `bcf_qc_to_parquet` — QC BCF WGS genotypes and convert to parquet

Runs genotype- and variant-level QC on a batch of BCF/VCF files and
consolidates the results into a long genotype and variant metadata parquet files.

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
  -iinput_file_list="/ukbgym_file_lists/ukbgym_gene_vcf_files.parquet" \
  -iaf_threshold=0.001 \
  -ifasta_ref="/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  -ifasta_ref_index="/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" \
  -isamples_file="/processed_data/phenotypes/eur_samples.txt" \
  --destination /processed_data/genotypes/
```

### `vep_loftee_parallel` — annotate variants with VEP + LOFTEE

Annotates the variant metadata with Ensembl VEP and the LOFTEE plugin. Input variants are split into
chunks that are annotated in parallel subjobs. Uses `variant_metdata.parquet` created in the `bcf_qc_to_parquet` stp as input. 


```bash
dx run vep_loftee_parallel \
  -ivariants_parquet="/processed_data/genotypes/variant_metadata.parquet" \
  -ivep_version="113" \
  -ichunk_size=500000 \
  -iuse_loftee=true \
  -igenes_to_keep_file="/ukbgym_file_lists/ukbbgym_gene_list.txt"
  --destination /processed_data/annotations/
```

---


## Prepare olink proteomics data 

### `prepare_regenie_olink_inputs` — format Olink data for REGENIE
<!-- TODO add the file list olink_covariate_fields.txt for olink and olink_covariate_field_mapping.parquet-->

Extract covariates and protein levels for Olink data and write it to the three tab-separated files REGENIE expects.
Can be run independently of all other steps. 
#### 1. Export the inputs from the dataset first with `table-exporter`:

<!-- TODO it seems like we need to run table exporter twice, once for covariates and once for proteomics levels?  -->

Extract covariates: 
```bash
dx run table-exporter \
  -idataset_or_cohort_or_dashboard=record-REDACTED \
  -ientity="participant" \
  -ifield_names_file_txt="/ukbgym_file_lists/olink_covariate_fields.txt" \
  -ioutput="olink_covariates_export" \
  -icoding_option="RAW" \
  -iheader_style="FIELD-NAME" \
  --destination="/processed_data/olink/raw/" \
  --priority normal \
  --instance-type="mem1_ssd1_v2_x16"
```

<!-- TODO check if this is correct -->
Extract OLINK protein levels: 
```bash
dx run table-exporter \
  -idataset_or_cohort_or_dashboard=record-REDACTED \
  -ientity="participant" \
  -ifield_names_file_txt="/ukbgym_file_lists/olink_fields.txt" \
  -ioutput="olink_export" \
  -icoding_option="RAW" \
  -iheader_style="FIELD-NAME" \
  --destination="/processed_data/olink/raw/" \
  --priority normal \
  --instance-type="mem1_ssd1_v2_x16"
```

#### 2. Run the preparation applet

```bash
dx run prepare_regenie_olink_inputs \
  -iraw_covariates="/processed_data/olink/raw/olink_covariates_export.csv" \
  -icovariate_mapping="/ukbgym_file_lists/olink_covariate_field_mapping.parquet" \
  -iproteomics_levels="/processed_data/olink/raw/olink_export.csv" \
  --destination /processed_data/olink/regenie/
```

### Run REGNIE step 1 on Olink proteomics data 
<!-- TODO explain how to run REGENIE Step 1  -->

Run REGENIE step 1 on these files to produce the `.prs` files needed by the next applet.



### `adjust_olink_ukbgym` — PRS-adjust Olink proteomics
<!-- TODO explain how do we get to /processed_data/olink/samples.txt?  I guess we actually don't need this but can just inner join with ukbgym samples? Or was there any other sample-level QC? -->

Normalises the Olink protein levels, regresses out the REGENIE polygenic risk scores per protein (OLS),
and run [Protrider](https://github.com/gtsitsiridis/protadjust) — an autoencoder — on the residuals to
flag outliers. 

Requries 

```bash
dx run adjust_olink_ukbgym \
  -iinput_csv="/processed_data/olink/raw/olink_export.csv" \
  -iregenie_step_1_prs_list="/regenie/step1_prs.list" \
  -isample_list="/ukbgym_file_lists/olink_samples.txt" \
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
