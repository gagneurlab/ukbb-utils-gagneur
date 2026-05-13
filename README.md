# ukbb-utils Gagneur

Pipeline for running [DeepRVAT](https://github.com/PMBio/deeprvat) on the UK Biobank Research Analysis Platform (RAP) using Whole Genome Sequencing (WGS) data.

## Overview

This repository contains DNAnexus applets and preprocessing scripts for the full DeepRVAT-WGS analysis workflow, from raw WGS BCF files through to association testing results.

```
WGS BCF files (RAP)
    └─[scripts/wgs_qc]──────────────── QC'd BCF (protein-coding genes ±10kb)
         └─[applets/vcf_qc_to_parquet] gt_long.parquet + variant_metadata.parquet
              └─[applets/vep_loftee_parallel]──── VEP-annotated variants
                   └─[applets/add_variant_annotations]── fully annotated variants
                                                              │
UKB Phenotypes (RAP dataset)                                  │
    └─[applets/extract_phenotypes_and_covariates]─────────────┤
         └─[applets/association_testing]────────────── association results

Olink Proteomics (optional)
    └─[applets/prepare_regenie_olink_inputs]
         └─[applets/adjust_olink_ukbgym]──── adjusted proteomics
```

## Repository Structure

```
deeprvat_wgs/
├── applets/          # DNAnexus applets (cloud-executable pipeline steps)
├── scripts/          # Preprocessing scripts and analysis notebooks
├── utils.py          # Shared DNAnexus API utilities
├── plotting_config.yaml
├── CITATION.cff
└── LICENSE
```

## Setup

### Pre-commit hooks

**Required:** install pre-commit hooks so that Jupyter notebook outputs are stripped and formatting is applied automatically on every commit.

```bash
pip install pre-commit
pre-commit install
```

### DNAnexus CLI

All applets run on the UK Biobank RAP. Install and configure the DNAnexus CLI:

```bash
pip install dxpy
dx login
```

See [scripts/rap_usage/README.md](scripts/rap_usage/README.md) for a full RAP quickstart guide covering SSH, VSCode, dxfuse, CONDA environments, and job cost monitoring.

## Applets

See [applets/README.md](applets/README.md) for descriptions and input/output specs for each applet.

| Applet | Purpose |
|--------|---------|
| [vcf_qc_to_parquet](applets/vcf_qc_to_parquet/) | QC BCF/VCF files and convert to parquet |
| [vep_loftee_parallel](applets/vep_loftee_parallel/) | VEP + LOFTEE variant annotation |
| [add_variant_annotations](applets/add_variant_annotations/) | Merge additional effect prediction scores |
| [extract_phenotypes_and_covariates](applets/extract_phenotypes_and_covariates/) | Extract UKB phenotypes and run REGENIE |
| [association_testing](applets/association_testing/) | Burden t-tests with covariate adjustment |
| [hpopt_training](applets/hpopt_training/) | DeepRVAT model training with hyperparameter optimisation |
| [prepare_regenie_olink_inputs](applets/prepare_regenie_olink_inputs/) | Format Olink proteomics for REGENIE |
| [adjust_olink_ukbgym](applets/adjust_olink_ukbgym/) | PRS-adjust Olink proteomics with Protrider |

## Scripts

See [scripts/README.md](scripts/README.md) for details on each preprocessing component.

## Citation

If you use this software, please cite it using the metadata in [CITATION.cff](CITATION.cff).
