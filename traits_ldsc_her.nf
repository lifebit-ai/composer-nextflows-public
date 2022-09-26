#!/usr/bin/env nextflow
/*
========================================================================================
                         bi-traits-nf
========================================================================================
 bi-traits-nf Analysis Pipeline.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*---------------------------------------
  Define and show help message if needed
-----------------------------------------*/

def helpMessage() {

    log.info"""
    
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --analysis_mode heritability --input_gwas_statistics GWAS123.vcf
    Essential parameters:
    --analysis_mode                  Type of analysis desired. Options are 'heritability' or 'genetic_correlation'.
    --input_gwas_statistics          Path to GWAS summary statistics file in GWAS VCF format.
    
    Optional parameters:
    --method                         Software used to perform the analysis. Default = LDSC. Currently available input options: LDSC, LDAK, GCTA_GREML.
    --other_gwas_statistics          Path to second set of GWAS summary statistics to be used for genetic correlation.                    
    --hapmap3_snplist                Path to SNP list from Hapmap needed for seleting SNPs considered for analysis
    --ld_scores_tar_bz2              Path to tar.bz2 files with precomputed LD scores. Alternatively, population can be specified via --pop parameter to use population-specific 1000Genomes LD scores. 
                                     If both --ld_scores_tar_bz2 and --pop are specified, LD scores provided via --ld_scores_tar_bz2 will be used.
    --pop                            Population (determines population-specific 1000Genomes LD scores used). Can be specified 
                                     instead of --ld_scores_tar_bz2 parameter. Default = EUR. Current available input options: EUR (European), EAS (East Asian), GBR (British).
                                     If both --ld_scores_tar_bz2 and --pop are specified, LD scores provided via --ld_scores_tar_bz2 will be used.
    --thin_ldak_tagging_file         Path to thin tagging model file used exclusively in the LDAK mode.
    --bld_ldak_tagging_file          Path to bld tagging model file used exclusively in the LDAK mode.
    --outdir                         Path to output directory
    --output_tag                     String containing output tag
    --gwas_sample_size               Number of samples in the input GWAS VCF (int)
                                     (Default: $params.gwas_sample_size)
    --other_gwas_sample_size      Number of samples in the external GWAS VCF (int)
                                     (Default: $params.other_gwas_sample_size)
    """.stripIndent()
}

// Show help message

if (params.help) {
    helpMessage()
    exit 0
}



/*---------------------------------------------------
  Define and show header with all params information 
-----------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Output dir']                     = params.outdir
summary['Launch dir']                     = workflow.launchDir
summary['Working dir']                    = workflow.workDir
summary['Script dir']                     = workflow.projectDir
summary['User']                           = workflow.userName

summary['analysis_mode']                  = params.analysis_mode
summary['input_gwas_statistics']          = params.input_gwas_statistics
summary['hapmap3_snplist']                = params.hapmap3_snplist
summary['ld_scores_tar_bz2']              = params.ld_scores_tar_bz2                        = params.method
summary['output_tag']                     = params.output_tag
summary['outdir']                         = params.outdir
summary['gwas_sample_size']               = params.gwas_sample_size

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*--------------------------------------------------
  LDSC - Genetic correlation and heritability
---------------------------------------------------*/
process prepare_vcf_files_ldsc {
    tag "preparation_files"
    label 'bcftools'
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(gwas_vcf)

    output:
    path("${params.output_tag}_transformed_gwas_stats.txt"), emit: ldsc_input

    script:

    """
    echo "CHR POS SNPID Allele1 Allele2 BETA SE p.value" > base.data.pre
    bcftools query -f'%CHROM %POS [%SNP] %REF %ALT [%BETA] [%SE] [%P]\n' $gwas_vcf >> base.data.pre
    # Generating the N column
    echo "N" > n_col.txt
    for i in \$(seq 2 `wc -l < base.data.pre`); do
      echo $params.gwas_sample_size >> n_col.txt
    done
    # Generating the imputationInfo column
    echo "imputationInfo" > info_col.txt
    for i in \$(seq 2 `wc -l < base.data.pre`); do
      echo "1" >> info_col.txt
    done
    paste -d " " base.data.pre n_col.txt info_col.txt > "${params.output_tag}_transformed_gwas_stats.txt"
    """
}

    
process reformat_for_ldsc {
      tag "reformat_for_ldsc"
      publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'
      stageInMode 'copy'

      input:
      file(ldsc_summary_stats)
      file(hapmap3_snplist)

      output:
      path("${params.output_tag}_ldsc.sumstats.gz"), emit: saige_ldsc

      script:

      """
      munge_sumstats.py --sumstats $ldsc_summary_stats \
                        --out "${params.output_tag}_ldsc" \
                        --merge-alleles $hapmap3_snplist \
                        --a1 Allele1 \
                        --chunksize ${params.munge_sumstats_chunksize} \
                        --a2 Allele2 \
                        --signed-sumstats BETA,0 \
                        --p p.value \
                        --snp SNPID \
                        --info imputationInfo
      """
}

process heritability_ldsc {
    tag "heritability"
    publishDir "${params.outdir}/heritability/", mode: 'copy'

    input:
    file(saige_output)
    file(ld_scores_tar_bz2)

    output:
    path("${params.output_tag}_h2.log"), emit: ldsc_report_input

    script:
    """
    tar -xvjf ${ld_scores_tar_bz2}
    ldsc.py \
      --h2 $saige_output \
      --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --out ${params.output_tag}_h2
    """
}

workflow lifebitai_traits_ldsc_her{
  take:
    ch_gwas_statistics
    ch_hapmap3_snplist
    ch_ld_scores_tar_bz2

  main:
    prepare_vcf_files_ldsc(ch_gwas_statistics)

    reformat_for_ldsc(prepare_vcf_files_ldsc.out.ldsc_input, 
                        ch_hapmap3_snplist)

    heritability_ldsc(reformat_for_ldsc.out.saige_ldsc,
                      ch_ld_scores_tar_bz2)
    
    ldsc_out =  heritability_ldsc.out.ldsc_report_input

  emit:
    ldsc_out

}

workflow {
  if (params.method == 'LDSC') {
    if (params.analysis_mode == 'heritability'){
      ch_gwas_statistics   =  Channel
                              .fromPath(params.input_gwas_statistics)
                              .ifEmpty { exit 1, "Cannot find GWAS stats file : ${params.input_gwas_statistics}" }
      
      ch_hapmap3_snplist   =  Channel
                              .fromPath(params.hapmap3_snplist)
                              .ifEmpty { exit 1, "Cannot find HapMap3 snplist file : ${params.hapmap3_snplist}" }

      ch_ld_scores_tar_bz2 = Channel
                              .fromPath(params.ld_scores_tar_bz2)
                              .ifEmpty { exit 1, "Cannot find LD scores file : ${params.ld_scores_tar_bz2}" }
      
      lifebitai_traits_ldsc_her(ch_gwas_statistics,
                                  ch_hapmap3_snplist,
                                  ch_ld_scores_tar_bz2)

    }
  }
}