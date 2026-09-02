# Extract Phenotypes and Covariates (DNAnexus Platform App)

Extracts phenotype and covariate fields directly from the UKB dataset. Optionally subsets to European ancestry and unrelated individuals, and apply statin correction to cholesterol/LDL.

If PLINK `.bed`/`.bim`/`.fam` files are supplied, the applet also runs REGENIE and returns polygenic
risk scores and PRS-corrected phenotypes; without them, only phenotypes and covariates are produced.


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


Alongside the phenotype parquets it writes two sample lists, one eid per line, no header:

| Output | Cohort | Feed it to |
|---|---|---|
| `eur_samples.txt` | after the ancestry filter, **before** relatedness pruning | `bcf_qc_to_parquet` → `samples_file` |
| `eur_unrelated_samples.txt` | after ancestry filtering **and** relatedness pruning | `avg_pheno_per_variant_traits` → `samples_txt` |

Both are snapshotted before the statin correction and quantile transform, and before the PRS
join, which drops participants without array genotype calls. With `subset_unrelated=False` the
two files are identical; with `subset_eur=False` the first is the full extracted cohort.


```bash
# upload the field lists once; reuse the file IDs on later runs
pheno_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_trait_fieldIDs.txt --brief)
covar_file=$(dx upload applets/extract_phenotypes_and_covariates/ukbbgym_covariate_fieldIDs.txt --brief)

dx run extract_phenotypes_and_covariates \
  -idataset_id="record-XXXXXXXX" \
  -ipheno_list_file="$pheno_file" \
  -iadditional_covars_file="$covar_file" \
  -ibed_file="/plink/ukb.bed" -ibim_file="/plink/ukb.bim" -ifam_file="/plink/ukb.fam" \
  -isubset_eur=true \
  -isubset_unrelated=true \
  -istatin_correction=true \
  --destination /processed_data/phenotypes/
```
