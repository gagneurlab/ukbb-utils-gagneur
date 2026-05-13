#!/bin/bash
set -e

OUTPUT_DIR="/regenie_step1_data"

echo "Submitting Phase 1: Merging Genotypes..."
# We use --brief to capture ONLY the job ID (e.g., job-Gxxxxx89...)
MERGE_JOB_ID=$(dx run app-swiss-army-knife \
  --instance-type="mem2_ssd1_v2_x16" \
  --destination="${OUTPUT_DIR}" \
  -icmd='
    set -e -x
    echo "1. Downloading raw chromosomes..."
    dx download "$DX_PROJECT_CONTEXT_ID:/Bulk/Genotype Results/Genotype calls/ukb22418_c*_b0_v2.*"
    
    echo "2. Merging..."
    ls ukb22418_c*.bed | sed "s/.bed//" > merge_list.txt
    FIRST=$(head -n 1 merge_list.txt)
    sed -i "1d" merge_list.txt
    
    plink --bfile "$FIRST" --merge-list merge_list.txt --make-bed --out ukb_cal_allChrs
    
    # Clean up so only the merged files are uploaded
    rm ukb22418_c* merge_list.txt
  ' --brief -y)

echo "Merge job submitted! ID: $MERGE_JOB_ID"
echo "Submitting Phase 2: QC and Pruning (Will wait for Merge to finish)..."

# Notice the --depends-on flag below!
dx run app-swiss-army-knife \
  --depends-on "$MERGE_JOB_ID" \
  --instance-type="mem2_ssd1_v2_x16" \
  --destination="${OUTPUT_DIR}" \
  -icmd='
    set -e -x
    PREFIX="ukb_cal_allChrs"
    
    echo "1. Downloading merged dataset from project..."
    # We download via dx command to avoid submission-time file validation errors
    dx download "$DX_PROJECT_CONTEXT_ID:'${OUTPUT_DIR}'/${PREFIX}.*"

    echo "2. Applying Strict QC..."
    plink2 --bfile $PREFIX --maf 0.01 --mac 100 --geno 0.1 --mind 0.1 --hwe 1e-15 \
           --write-snplist --write-samples --no-id-header --out step1_filtered

    echo "3. LD Pruning..."
    plink2 --bfile $PREFIX --keep step1_filtered.id --extract step1_filtered.snplist \
           --indep-pairwise 1000 100 0.9 --out step1_pruned

    echo "4. Generating Final Dataset..."
    plink2 --bfile $PREFIX --keep step1_filtered.id --extract step1_pruned.prune.in \
           --make-bed --out regenie_step1_ready

    # Clean up so only the finalized REGENIE files are uploaded
    rm ${PREFIX}* step1_filtered* step1_pruned*
  ' -y

echo "Pipeline submitted successfully! You can close your terminal."