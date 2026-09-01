<!-- dx-header -->
# Prepare regenie Olink inputs (DNAnexus Platform App)

Formats raw UKB covariates and Olink protein levels into the three tab-separated files REGENIE expects.

Maps raw UKB field IDs to readable labels via `olink_covariate_field_mapping.parquet`, one-hot encodes
categorical covariates, adds derived features (age², age×sex, …), then joins the covariates to the
protein levels on sample ID and writes out covariates, phenotypes and a sample list.

Both CSV inputs must be exported from the dataset with `table-exporter` before running this applet.
Run REGENIE step 1 on the outputs to produce the `.prs` files needed by `adjust_olink_ukbgym`.
<!-- /dx-header -->

Inputs, outputs and example commands: see the [repository README](../../README.md).
