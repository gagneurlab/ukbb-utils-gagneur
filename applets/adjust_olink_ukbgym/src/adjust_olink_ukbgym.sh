#!/bin/bash
# protadjust 0.0.1

set -euxo pipefail

main() {
    # install uv
    curl -LsSf https://astral.sh/uv/install.sh | sh

    # install protadjust
    PROTADJUST_IMAGE="ghcr.io/gtsitsiridis/protadjust:latest"
    echo "Pulling Docker image $PROTADJUST_IMAGE ..."
    docker pull "$PROTADJUST_IMAGE"
    protadjust() {
        docker run --rm -v "$PWD:/data" "$PROTADJUST_IMAGE" "$@"
    }

    # input_csv="file-REDACTED"
    # sample_list="file-REDACTED"
    # regenie_step_1_prs=($(dx ls project-REDACTED:/processed_data/olink/regenie/step_1/*.prs --brief))
    # regenie_step_1_prs_list=$(dx ls project-REDACTED:/processed_data/olink/regenie/step_1/*_prs.list --brief)

    # --- download inputs ---
    mkdir -p csv_files
    for file_id in "${input_csv[@]}"; do
        dx download "$file_id" --output csv_files/ &
    done
    wait
    first=1
    for f in csv_files/*.csv; do
        if [ $first -eq 1 ]; then
            cat "$f"
            first=0
        else
            tail -n +2 "$f"
        fi
    done > input.csv
    dx download "$sample_list" -o sample_list.txt
    dx download "$regenie_step_1_prs_list" -o prs.list
    mkdir -p prs_files
    for file_id in "${regenie_step_1_prs[@]}"; do
        dx download "$file_id" --output prs_files/ &
    done
    wait
    wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/olink_assay.dat

    # Step 1: preprocess olink
    uv run preprocess_olink.py --olink-path input.csv --helper-assay-path olink_assay.dat --samples-path sample_list.txt --output-dir .

    # Step 2: aggregate PRS files into parquet
    uv run aggregate_prs.py --prs-list prs.list --prs-dir prs_files/ --output prs.parquet

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
