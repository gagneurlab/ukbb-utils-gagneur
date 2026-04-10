<!-- dx-header -->
# Adjust Olink proteomics (DNAnexus Platform App)

Preprocess Olink proteomics data, adjust for polygenic risk scores (PRS) via per-protein regression, and detect outliers using the Protrider autoencoder.
<!-- /dx-header -->

## Inputs

| Name | Type | Required | Description |
|---|---|---|---|
| `input_csv` | file | Yes | Olink CSV file exported using TableExporter |
| `regenie_step_1_prs` | array:file | No | `.prs` files from regenie step 1 output |
| `regenie_step_1_prs_list` | file | No | `.list` file mapping gene names to `.prs` filenames (space-separated, no header) |
| `sample_list` | file | No | Line-separated text file of sample IDs to include |

## Output

| Name | Type | Description |
|---|---|---|
| `output_parquet` | file | Protrider-adjusted protein abundance values (`adjusted_proteomics.parquet`) |

## Prerequisites

`input_csv` must be exported from the UKB RAP dataset before running this applet. Use the `table-exporter` app:

```bash
dx run table-exporter \
  -idataset_or_cohort_or_dashboard=record-REDACTED \
  -ientity="olink_instance_0" \
  -ioutput="olink_export" \
  -icoding_option="RAW" \
  -iheader_style="FIELD-NAME" \
  --destination="/processed_data/olink/raw/" \
  --priority normal \
  --instance-type="mem1_ssd1_v2_x16"
```

This produces the CSV at `/processed_data/olink/raw/olink_export.csv` on the platform. Pass that file as `input_csv`.

## Pipeline steps

1. **Preprocess Olink** — convert CSV to parquet, optionally filter to `sample_list`, normalise using the Olink assay helper file
2. **Aggregate PRS** — download all `.prs` files and combine into a single `prs.parquet` (samples × proteins)
3. **PRS regression** — per-protein OLS regression against PRS scores; output is residuals parquet
4. **Protrider** — autoencoder-based outlier detection on the PRS residuals

## Instance type

Default instance: `mem1_ssd1_v2_x16` (16 cores, ~120 GB RAM).
