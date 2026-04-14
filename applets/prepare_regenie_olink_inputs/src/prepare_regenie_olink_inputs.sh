#!/bin/bash
# prepare_regenie_olink_inputs 0.0.1

set -euxo pipefail

main() {
    # install uv
    curl -LsSf https://astral.sh/uv/install.sh | sh
    source $HOME/.local/bin/env
    
    dx-download-all-inputs --parallel

    uv run preprocess_covariates.py \
        --covariates-path "$raw_covariates_path" \
        --mapping-path "$covariate_mapping_path" \
        --output-dir .

    uv run prepare_regenie_olink_inputs.py \
        --covariates-parquet covariates_preprocessed.parquet \
        --levels-csv "$proteomics_levels_path" \
        --output-dir .

    regenie_covariates=$(dx upload regenie_proteomics_covariates.txt --brief)
    dx-jobutil-add-output regenie_covariates "$regenie_covariates" --class=file

    regenie_levels=$(dx upload regenie_proteomics_levels.txt --brief)
    dx-jobutil-add-output regenie_levels "$regenie_levels" --class=file

    regenie_samples=$(dx upload regenie_proteomics_samples.txt --brief)
    dx-jobutil-add-output regenie_samples "$regenie_samples" --class=file
}
