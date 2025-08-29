# deeprvat_wgs

Repository having code for running DeepRVAT2 on the UK Biobank Research Analysis Platform (RAP) using Whole Genome Sequencing (WGS) data.

## Usage !! Setup pre-commit!!

**!! Setup pre-commit!!** so formatting and removal of jupyter notebook outputs is automatically done.
To do this, just install `pre-commit` into your environment and run `pre-commit install` in your repository after you have cloned it for the first time. Then, pre-commit will always run when you commit.

## [`wgs_qc/`](./wgs_qc/)

Provides code for extracting the regions covered by each WGS file on the RAP and then runs QC on each file.

- The file providing the **extracted regions** for each bcf file is here: `/s/project/deeprvat/wgs_preprocessing_rap/qc/bcf_region_files/all_region_files.parquet` - The **qced bcf files** (covering all protein coding genes +- 10kb) can be found here (RAP): `/ukb-gagneur/processed_data/wgs/bcf_files_qced/batch1/`.

QC applied on each file (as in [Hawkes et al.](https://www.nature.com/articles/s41588-025-02095-4#Sec11))

```
bcftools +setGT --output-type u -- -t q -i "FMT/GQ<=10 | smpl_sum(FMT/LAD)<8" -n . |
bcftools filter --output-type u -e "F_MISSING > 0.1" |
```

## [`vcf_2_zarr/`](./vcf_2_zarr/)

Includes DNAnexus applets ([`vcf2zarr_applet`](./vcf_2_zarr/vcf2zarr_applet/Readme.md) and [`vcf2zarr_encode_applet`](./vcf_2_zarr/vcf2zarr_encode_applet/Readme.md)) designed to convert VCF files into the Zarr format using the `vcf2zarr` tool within the RAP environment. These applets facilitate the partitioning and encoding steps and adapt them for the RAP.

## [`rap_usage/`](./rap_usage/README.md)

Provides a guide and useful commands for interacting with the RAP environment, including setting up environments, accessing data, running jobs, and using tools like VSCode and dxfuse. See the [RAP usage README](./rap_usage/README.md) for a quickstart.

## [`annotation/`](./annotation/)

Contains scripts and notebooks for preparing and merging various variant annotations required by DeepRVAT. See the [annotation README](./annotation/README.md) for details on the specific annotation sources and merging process.
