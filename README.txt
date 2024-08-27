# This is a repository with the combined ancestry processing steps for AGD data as a single WDL 

This WDL will run ancestry estimation (PCA & SCOPE) on AGD data. It runs in a series of 3 optional steps, which can be opted in to or out of using Boolean variables: 

1. Combine AGD data with Spike-In data (Boolean external_spike_in)
2. Run PCA projection from the GBMI group on the combined data OR on the AGD data alone, if the Spike-In data is not provided (Boolean run_pca)
3. Run unsupervised SCOPE on the combined data OR on the AGD data alone, if the Spike-In data is not provided (Boolean run_scope)
4. Run supervised SCOPE on the combined data OR on the AGD data alone, if the Spike-In data is not provided. (Boolean scope_supervised)
    4b. For step 4, the user can provide a file with within superpopulation frequencies, or alternatively the frequencies can be calculated from provided plink files, such as the Spike-In population.

# Description of input data 

Start with a set of AGD chromosomes in a data table, such that the following columns are arrays of strings or of file paths: chromosomes, pgen, psam and pvar files. 
Include only the autosomes. 

## Required input data 

- id_map_file (File): AGD ID map file: "gs://working-set-redeposit/ica-agd/cohort_001/20240303_agd35k_ica_primary_eligible.txt"
- Chromosomes (Array[String]): chromosomes to process: this.agd35k_bed_alls.chromosome
- source_pgen_files (Array[File]): AGD pgen files: this.agd35k_bed_alls.pgen_pgen
- source_psam_files s  (Array[File]): AGD psam files: this.agd35k_bed_alls.pgen_psam
- source_pvar_files  (Array[File]): AGD pvar files: this.agd35k_bed_alls.pgen_pvar

- target_prefix (String): 20240827_AGD35K_ancestry - prefix for all the output files 

- target_gcp_folder (String):  GCP folder to which the output files will be copied: "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/AGD_ancestry_pipeline/"

## Required input choices that have defaults: 

- external_spike_in (Boolean): true
- run_pca (Boolean): true
- run_scope (Boolean): true
- scope_supervised (Boolean): true


## AGD combine with spike-in 

These are optional inputs for the Spike-In population, but they are REQUIRED if running with external_spike_in=true. Here, the examples are for the TGP population (thousand genomes). 

Because pmerge does not yet work to merge pgen files, we have to use the Plink1 format, so if the user wants to use the Spike-In, they must provide the Plink1 version files (bed/bim/fam). This is consequently the most expensive part of the pipeline, since it requires large disk sizes and memory 

- spike_in_relatives_exclude (File): a list of IDs to exclude so that close relatives aren't included in the SCOPE estimation with these samples: "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/TGP_helper_files/deg2_hg38.king.cutoff.out.id"
- spike_in_pgen_file: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.pgen"
- spike_in_psam_file: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.psam"
- spike_in_pvar_file: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.pvar"
- Array[File]? source_bed_files
- Array[File]? source_bim_files
- Array[File]? source_fam_files

## PCA Inputs 

These are optional inputs for the PCA, but they are REQUIRED if running with run_pca=true. Here, the PCA is run using the variants from the TGP + HGDP populations as provided by GBMI. 
Since the AGD data is GRCh38, the PCA variants, loadings, and allele frequencies should also be in GRCh38.
The data for the variants, loadings, and allele frequencies we used can be accessed here: https://github.com/globalbiobankmeta/pca_projection/blob/master/doc_templates/prerequisites.md

- pca_variants_extract_file "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/pca_ref_files/variants.extract"
- pca_loadings_file (File): Loadings for the PCs as provided by GBMI. "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/pca_ref_files/hgdp_tgp_pca_gbmi_snps_loadings.GRCh38.plink.tsv"
- pca_af_file:  "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/pca_ref_files/hgdp_tgp_pca_gbmi_snps_loadings.GRCh38.plink.afreq"

## SCOPE inputs - Unsupervised 

These are optional inputs for running SCOPE, but they are REQUIRED if running with run_scope=true. SCOPE can always be run in the unsupervised mode. Optionally, SCOPE can also be run in a supervised mode, if reference population allele frequencies are provided. 

- Long Range LD files: this is optional in the WDL inputs because you may not run SCOPE but if you DO RUN SCOPE, this is required. "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/scope_ancestry_estimation/grch38_ld_long_range.txt"
- OPTIONAL:  scope_plink2_maf_filter: this is optional because the default will be used (maf 0.01): "--maf 0.01"
- OPTIONAL: scope_plink2_LD_filter_option: this is optional because the default will be used ("--indep-pairwise 50000 80 0.1"): "--indep-pairwise 50000 80 0.1"
- OPTIONAL: K (Int) - number of ancestry groups to use in SCOPE estimation. Default is 4, we used 5 to match the number in the reference supervised approach for TGP reference population (5 superpopulation continental ancestry) 
- OPTIONAL: seed (Int) - random seed for SCOPE.

## SCOPE inputs - Supervised

These are optional inputs for running SCOPE in the supervised mode, but they are REQUIRED if running with scope_supervised=true. Here, we used the TGP (same as spike-in) to calculate reference allele frequencies within superpopulations. 
The user can provide a file with within superpopulation frequencies, or alternatively the frequencies can be calculated from provided plink files, such as the Spike-In population.
To successfully run SCOPE in the supervised mode, the user must provide either the within superpopulation frequencies or the Spike-In population plnk files.
Here, we used the TGP (same as spike-in) to calculate reference allele frequencies within superpopulations.

- OPTIONAL: supervised_scope_reference_freq: within superpopulation frequency file. This is optional because the user could chose to calculate the superpopulation allelic frequencies within the WDL using provided plink files instead. 
- OPTIONAL: supervised_scope_reference_pgen_file: This is optional because the user could choose to provide precalculated within superpopulation allele frequencies instead: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.pgen"
- OPTIONAL: supervised_scope_reference_pvar_file: This is optional because the user could choose to provide precalculated within superpopulation allele frequencies instead: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.pvar"
- OPTIONAL: supervised_scope_reference_psam_file: This is optional because the user could choose to provide precalculated within superpopulation allele frequencies instead: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.psam"
- OPTIONAL: supervised_scope_reference_superpop_file: file with superpopulation groups for the supervised reference population: This is optional because the user could choose to provide precalculated within superpopulation allele frequencies instead: "gs://fc-2ee2ca2a-a140-48a1-b793-e27badb7945d/high_coverage_3202_samples/all_hg38_ns.superpopulations.txt"
- OPTIONAL (task option, not workflow option): plink2_maf_filter: this is optional because there is a default of 0.001. This is the maf filter for which SNPs will get allelic frequencies calculated within super population  ("--maf 0.001"). We used --maf 0.05

## Inputs that can be changed throughout 

- docker: the docker image to use for the tasks. Do not alter if you are not sure what you are doing, these have the correct docker image for each task as the default 
- memory_gb: the memory to use for each task. Alter if needed 
- out_prefix for some tasks: not recommended to change 