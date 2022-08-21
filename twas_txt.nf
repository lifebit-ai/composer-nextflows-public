#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --gwas_summary_statistics gwas_summary.tsv.gz [Options]
    
    Inputs Options:
    --gwas_summary_statistics        Path to input GWAS summary statistics file.
    --ld_reference_panel             Path to LD reference panel.
    --eqtl_weights                   Path to eQTL weights.
    
    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}



process ptwas_scan {
    tag "ptwas_scan"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(vcf_sumstats)
    file(ld_reference_panel)
    file(eqtl_weights)
    
    output:
    tuple path("*stratified_out.txt"), path("*summary_out.txt"), emit: gambit_output

    script:
    """
    tar xvzf ${ld_reference_panel}
    tabix -p vcf -f ${eqtl_weights}
    tabix -p vcf -f ${vcf_sumstats}
    ${params.gambit_exec_path} --gwas ${vcf_sumstats} --betas ${eqtl_weights} --ldref G1K_EUR_3V5/chr*.vcf.gz --ldref-only 
    """
  }

workflow lifebit_twas_txt{

    take:
        ch_ld_reference
        ch_eqtl_weights
        ch_gwas_sumstats

    main:
        if (params.help) {
            helpMessage()
            exit 0
        }

        ptwas_scan(ch_gwas_sumstats,
                    ch_ld_reference,
                    ch_eqtl_weights)

    emit:
        ptwas_scan.out.gambit_output
}
