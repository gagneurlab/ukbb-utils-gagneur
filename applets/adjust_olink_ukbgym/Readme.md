<!-- dx-header -->
# Adjust Olink proteomics (DNAnexus Platform App)

Normalises Olink protein levels, regresses out the REGENIE polygenic risk scores per protein (OLS), and
runs Protrider — an autoencoder — on the residuals to flag outliers.

Steps: preprocess the Olink CSV (optionally subset to `sample_list`) → aggregate the `.prs` files into
one samples × proteins table → per-protein OLS regression against the PRS → Protrider on the residuals.

The Olink CSV must be exported from the dataset with `table-exporter` (entity `olink_instance_0`) before
running this applet.
<!-- /dx-header -->

Inputs, outputs and example commands: see the [repository README](../../README.md).

## Developer notes

The applet is a thin shell wrapper around the
[`ghcr.io/gtsitsiridis/protadjust`](https://github.com/gtsitsiridis/protadjust) Docker image. Each run
pulls the image, downloads the inputs, runs `protadjust` with the working directory mounted at `/data`,
and uploads `output/adjusted_proteomics.parquet`. The image is pinned to `:latest` via
`PROTADJUST_IMAGE` in [`src/adjust_olink_ukbgym.sh`](src/adjust_olink_ukbgym.sh), so a new image release
does not require rebuilding the applet.

Protrider is memory-hungry. The default instance is `mem1_ssd1_v2_x16`; override it at launch with
`--instance-type mem2_ssd1_v2_x32` if a run runs out of memory.
