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

summary['Max Resources']                  = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']                     = params.outdir
summary['Launch dir']                     = workflow.launchDir
summary['Working dir']                    = workflow.workDir
summary['Script dir']                     = workflow.projectDir
summary['User']                           = workflow.userName

summary['analysis_mode']                  = params.analysis_mode
summary['input_gwas_statistics']          = params.input_gwas_statistics
summary['other_gwas_statistics']       = params.other_gwas_statistics
summary['hapmap3_snplist']                = params.hapmap3_snplist
summary['ld_scores_tar_bz2']              = params.ld_scores_tar_bz2
summary['pop']                            = params.pop
summary['method']                         = params.method
summary['input_ldak_statistics']          = params.input_ldak_statistics
summary['thin_ldak_tagging_file']         = params.thin_ldak_tagging_file
summary['bld_ldak_tagging_file']          = params.bld_ldak_tagging_file
summary['output_tag']                     = params.output_tag
summary['outdir']                         = params.outdir
summary['gwas_sample_size']               = params.gwas_sample_size
summary['other_gwas_sample_size']      = params.other_gwas_sample_size

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

  
process gcta_calculate_grm {
    label 'gcta'
    publishDir "${params.outdir}/GCTA/GRM/", mode: 'copy'

    input:
    tuple val(plink_prefix), file(bed), file(bim), file(fam)

    output:
    path("GRM*"), emit: gcta_grm

    script:
    """
    gcta_v1.94.0Beta_linux_kernel_2_x86_64_static \
                            --bfile $plink_prefix \
                            --autosome \
                            --maf ${params.maf_cutoff} \
                            --grm-cutoff ${params.grm_cutoff} \
                            --make-grm \
                            --out GRM 
    """


}

process gcta_greml_h2 {
    label 'gcta'
    publishDir "${params.outdir}/GCTA/heritability", mode: 'copy'

    input:
    file(pheno_file)
    file("*")

    output:
    path("GCTA_GREML_snp_h2*"), emit: gcta_greml

    script:
    """
    gcta_v1.94.0Beta_linux_kernel_2_x86_64_static \
            --grm GRM \
            --reml \
            --pheno $pheno_file \
            --out GCTA_GREML_snp_h2
            """
}

process gcta_greml_gc {
    label 'gcta'
    publishDir "${params.outdir}/GCTA/genetic_correlation", mode: 'copy'

    input:
    file(pheno_file)
    file("*")

    output:
    path("GCTA_GREML_gc*"), emit: gcta_greml

    script:
    """
    gcta_v1.94.0Beta_linux_kernel_2_x86_64_static \
    --grm GRM \
    --reml-bivar \
    --pheno $pheno_file \
    --out GCTA_GREML_gc
    """
}


workflow lifebitai_traits_gcta_greml{
  take:
    ch_plink_gcta_input
    ch_pheno_file

  main:
    gcta_calculate_grm(ch_plink_gcta_input)

    if (params.analysis_mode == 'heritability') {
      gcta_greml_h2(ch_pheno_file, gcta_calculate_grm.out.gcta_grm)
      gcta_greml_out = gcta_greml_h2.out.gcta_greml
    } else if (params.analysis_mode == 'genetic_correlation') {
      gcta_greml_gc(ch_pheno_file, gcta_calculate_grm.out.gcta_grm)
      gcta_greml_out = gcta_greml_gc.out.gcta_greml
    }

  emit:
    gcta_greml_out
}

workflow {
  if (params.method == 'GCTA_GREML') {
    ch_plink_gcta_input = Channel
      .fromFilePairs("${params.gcta_plink_input}", size:3, flat : true)
      .ifEmpty { exit 1, "Input files in plink format not found: ${params.gcta_plink_input}." }
    ch_pheno_file  = Channel
      .fromPath("${params.pheno_file}")
      .ifEmpty { exit 1, "Phenotype file not found: ${params.pheno_file}."}

   lifebitai_traits_gcta_greml(ch_plink_gcta_input, ch_pheno_file)
  }
}