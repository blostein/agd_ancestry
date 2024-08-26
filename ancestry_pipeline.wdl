version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/genotype/Utils.wdl" as http_GenotypeUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/agd/AgdUtils.wdl" as http_AgdUtils


# WORKFLOW

workflow agd_ancestry_workflow{
    input{
        # workflow choices 

        Boolean external_spike_in = true 
        Boolean scope_supervised = true
        Boolean run_pca = true
        Boolean run_scope= true

        # required inputs original data as array of chromosomes 

        Array[File] source_pgen_files
        Array[File] source_pvar_files
        Array[File] source_psam_files
        Array[String] chromosomes
        String target_prefix

        File id_map_file

        # optional inputs for PCA - required if running PCA 

        File? pca_variants_extract_file

        # optional inputs for spike in data - required if merging spike in data 

        File? spike_in_pgen_file
        File? spike_in_pvar_file
        File? spike_in_psam_file
        File? spike_in_relatives_exclude 
    
        # optional inputs for supervised scope    - required if running supervised scope 

        File? supervised_scope_reference_pgen_file
        File? supervised_scope_reference_pvar_file
        File? supervsied_scope_reference_psam_file
        File? supervised_scope_reference_superpop_file   
        File? supervised_scope_reference_relatives_exclude
        File? supervised_scope_reference_freq

        #optional inputs for unsupervised scope - required if running unsupervised scope
        String? scope_plink2_maf_filter = "--maf 0.01"

        String? scope_plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
        File scope_long_range_ld_file
        Int K = 4
        Int seed = 1234


    }

    # If the user chose to use the supervised scope and there is no precalculated reference allele frequency provided, allele frequency can be calculated from provided plink files instead
    if(scope_supervised && !defined(supervised_scope_reference_freq)){
        call CalculateFreq{
        input: 
            pgen_file = supervised_scope_reference_pgen_file,
            pvar_file = supervised_scope_reference_pvar_file,
            psam_file = supervised_scope_reference_psam_file,
            superpop_file = supervised_scope_reference_superpop_file,
            relatives_exclude = supervised_scope_reference_relatives_exclude
        }
        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile{
                input:
                    source_file = CalculateFreq.freq_file,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }
    }

    # If the user chose to use an external spike in, then first the spike-in data must be merged with the original data, and then the pipeline can proceed 
    if(external_spike_in){
         scatter (idx in range(length(chromosomes))) {
            String chromosome = chromosomes[idx]
            File pgen_file = source_pgen_files[idx]
            File pvar_file = source_pvar_files[idx]
            File psam_file = source_psam_files[idx]

            call SubsetChromosomeTGP{
                input: 
                    pgen_file = spike_in_pgen_file,
                    pvar_file = spike_in_pvar_file,
                    psam_file = spike_in_psam_file,
                    chromosome = chromosome,
                    relatives_exclude = spike_in_relatives_exclude
            }

            call ConvertPgenToBed {
                input: 
                    pgen = pgen_file,
                    pvar = pvar_file,
                    psam = psam_file
            }

            call Merge1000genomesAGD {
                input:
                    agd_bed_file = PGenToBed.out_bed,
                    agd_bim_file = PGenToBed.out_bim,
                    agd_fam_file = PGenToBed.out_fam,
                    TGP_bed_file = SubsetChromosomeTGP.out_bed_file,
                    TGP_bim_file = SubsetChromosomeTGP.out_bim_file,
                    TGP_fam_file = SubsetChromosomeTGP.out_fam_file
            }
        }
    }

    #now we can proceed with the pipeline, using select_first to get the first non empty array value to select either the merged or original input files 

    Array[File] my_pgen_files=select_first([Merge1000genomesAGD.out_pgen_file, source_pgen_files])
    Array[File] my_pvar_files=select_first([Merge1000genomesAGD.out_pvar_file, source_pvar_files])
    Array[File] my_psam_files=select_first([Merge1000genomesAGD.out_psam_file, source_psam_files])

}

