<!-- dx-header -->
# Prepare regenie Olink inputs (DNAnexus Platform App)

Preprocesses raw UKB covariates and Olink proteomics levels into regenie-ready text files.
<!-- /dx-header -->

## Prerequisites

The raw covariates parquet must be exported from the UKB dataset using table-exporter with the bundled field list:

```bash
fields_file=$(dx upload olink_covariate_fields.txt --brief)

dx run table-exporter \
  -idataset_or_cohort_or_dashboard=record-REDACTED \
  -ientity="participant" \
  -ifield_names_file_txt="$fields_file" \
  -ioutput="olink_covariates_export" \
  -icoding_option="RAW" \
  -iheader_style="FIELD-NAME" \
  --destination="/processed_data/olink/raw/" \
  --priority normal \
  --instance-type="mem1_ssd1_v2_x16"
```

The field list (`olink_covariate_fields.txt`) is derived from `olink_covariate_field_mapping.parquet`.

## Inputs

| Name | Type | Required | Description |
|---|---|---|---|
| `raw_covariates` | file | Yes | Raw covariates CSV from table-exporter (samples x UKB fields) |
| `covariate_mapping` | file | Yes | Parquet mapping UKB field IDs to field names and titles (`olink_covariate_field_mapping.parquet`) |
| `proteomics_levels` | file | Yes | Proteomics levels CSV exported from table-exporter |

## Outputs

| Name | Type | Description |
|---|---|---|
| `regenie_covariates` | file | Tab-separated covariates file for regenie (FID, IID, covariates) |
| `regenie_levels` | file | Tab-separated phenotype file for regenie (FID, IID, proteins) |
| `regenie_samples` | file | Tab-separated sample list for regenie (FID, IID) |

## Pipeline steps

1. **Preprocess covariates** — map raw UKB fields to covariate labels, one-hot encode categoricals, compute derived features (age², ageXsex, etc.)
2. **Prepare regenie inputs** — join covariates with proteomics levels on sample ID, format FID/IID columns, write regenie text files

## Instance type

Default instance: `mem1_ssd1_v2_x2`. Scale up if the covariate matrix is large.
