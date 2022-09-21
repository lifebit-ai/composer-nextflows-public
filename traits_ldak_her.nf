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
summary['input_gwas_statistics']          = params.input_gwas_statistics
summary['bld_ldak_tagging_file']          = params.bld_ldak_tagging_file
summary['output_tag']                     = params.output_tag
summary['outdir']                         = params.outdir
summary['gwas_sample_size']               = params.gwas_sample_size

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

process prepare_vcf_ldak_summary {
  tag "preparation_files"
  label 'bcftools'
  publishDir "${params.outdir}/LDAK/vcf_input/", mode: 'copy'

  input:
  path(gwas_vcf) 

  output:
  path("ldak_transformed_input_sumstats.txt"), emit: ldak_transformed_input

  script:

  """
  echo "Predictor A1 A2 Direction P" > base.data.pre
  bcftools query -f'%CHROM:%POS %REF %ALT [%BETA] [%P]\n' $gwas_vcf >> base.data.pre
  # Generating the N column
  echo "n" > n_col.txt
  for i in \$(seq 2 `wc -l < base.data.pre`); do
    echo $params.gwas_sample_size >> n_col.txt
  done
  paste -d " " base.data.pre n_col.txt > "n_base_data.txt"

  # remove duplicate rows
  cat -n n_base_data.txt | sort -uk2 | sort -nk1 | cut -f2- > ldak_transformed_input_sumstats.txt
  """
}

process heritability_ldak {
  tag "heritability"
  label 'ldak'
  publishDir "${params.outdir}/LDAK/heritability/", mode: 'copy'

  input:
  file(ldak_summary_input)
  file(ldak_tagging)

  output:
  path("${params.output_tag}_snpher*"), emit: ldak_heritability_output

  script:
  """
  gunzip -f ${ldak_tagging}
  ldak5.2.linux \
    --sum-hers ${params.output_tag}_snpher \
    --summary ${ldak_summary_input} \
    --tagfile ${ldak_tagging.baseName} \
    --check-sums NO

  """
}

workflow lifebitai_traits_ldak_her{
  take:
    ch_gwas_statistics
    ch_bld_ldak_tagging_file

  main:
  prepare_vcf_ldak_summary(ch_gwas_statistics)
  
  heritability_ldak(prepare_vcf_ldak_summary.out.ldak_transformed_input,
                    ch_bld_ldak_tagging_file)

  ldak_out = heritability_ldak.out.ldak_heritability_output

  emit:
    ldak_out
}

workflow {
  if (params.method == 'LDAK') {
    /*--------------------------------------------------
  Channel preparation
---------------------------------------------------*/
    if (params.pop && params.ld_scores_tar_bz2){
      log.warn "Both LD scores and pre-computed 1000Genomes population-specific LD scores are provided. Custom LD scores are used."
    }
    if (params.analysis_mode == 'heritability'){
      ch_gwas_statistics        = Channel.
                                  fromPath(params.input_gwas_statistics)
                                  .ifEmpty { exit 1, "Cannot find GWAS stats file : ${params.input_gwas_statistics}" }

      ch_bld_ldak_tagging_file  = Channel.
                                  fromPath(params.bld_ldak_tagging_file)
                                  .ifEmpty { exit 1, "Cannot find bld LDAK tagging file : ${params.bld_ldak_tagging_file}" }

      lifebitai_traits_ldak_her(ch_gwas_statistics,
                            ch_bld_ldak_tagging_file)
    }
  }
}