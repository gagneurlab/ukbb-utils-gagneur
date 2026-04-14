#!/bin/bash
# protadjust 0.0.1

set -euxo pipefail

main() {
    # install uv
    curl -LsSf https://astral.sh/uv/install.sh | sh
    source $HOME/.local/bin/env

    # install protadjust
    PROTADJUST_IMAGE="ghcr.io/gtsitsiridis/protadjust:latest"
    echo "Pulling Docker image $PROTADJUST_IMAGE ..."
    docker pull "$PROTADJUST_IMAGE"
    protadjust() {
        docker run --rm -v "$PWD:/data" "$PROTADJUST_IMAGE" "$@"
    }

    # --- download inputs ---
    dx-download-all-inputs --parallel

    # concatenate input CSVs (array:file) into a single input.csv
    first=1
    for f in "${input_csv_path[@]}"; do
        if [ $first -eq 1 ]; then
            cat "$f"
            first=0
        else
            tail -n +2 "$f"
        fi
    done > input.csv

    wget -nd biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/olink_assay.dat

    # Step 1: preprocess olink
    uv run preprocess_olink.py \
        --olink-path input.csv \
        --helper-assay-path olink_assay.dat \
        --samples-path "$sample_list_path" \
        --output-dir .

    # Flatten array:file PRS inputs into a single directory (dx-download-all-inputs
    # places each file under a numbered subdirectory, e.g. .../regenie_step_1_prs/0/file.prs)
    mkdir -p prs_files
    for f in "${regenie_step_1_prs_path[@]}"; do
        cp "$f" prs_files/
    done

    # Step 2: aggregate PRS files into parquet
    uv run aggregate_prs.py \
        --prs-list "$regenie_step_1_prs_list_path" \
        --prs-dir prs_files \
        --output prs.parquet

    # Step 3: adjust proteomics data for PRS
    protadjust -v \
        protein-regression /data/olink_preprocessed.parquet /data/prs_residuals \
        --protein-covariate-path /data/prs.parquet --index-col sample --protein-covariate-index-col sample

    # Step 4: adjust PRS residuals using protrider
    protadjust -v protrider /data/prs_residuals/adjusted_proteomics.parquet /data/protrider \
        --pval-dist t --max-nas 0.3 --n-layers 1 --n-epochs 1000 --lr 0.001 --find-q-method OHT --index-col sample

    # --- upload output ---
    output_parquet=$(dx upload protrider/adjusted_proteomics.parquet --brief)
    dx-jobutil-add-output output_parquet "$output_parquet" --class=file
}
