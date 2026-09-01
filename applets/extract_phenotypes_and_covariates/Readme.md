<!-- dx-header -->
# Extract Phenotypes and Covariates (DNAnexus Platform App)

Extracts phenotype and covariate fields directly from the UKB dataset, batching requests to stay
under the API limits. Optionally subsets to European ancestry and unrelated individuals, and applies
statin correction to cholesterol/LDL.

If PLINK `.bed`/`.bim`/`.fam` files are supplied, the applet also runs REGENIE and returns polygenic
risk scores and PRS-corrected phenotypes; without them, only phenotypes and covariates are produced.

Field lists for the UKBGym analysis are kept in this directory alongside the applet source:
`ukbbgym_trait_fieldIDs.txt` and `ukbbgym_covariate_fieldIDs.txt`. They are not shipped with the built
applet — upload them to the project with `dx upload` and pass the file IDs as `pheno_list_file` and
`additional_covars_file`.
<!-- /dx-header -->

Inputs, outputs and example commands: see the [repository README](../../README.md).
