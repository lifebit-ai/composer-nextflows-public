#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process transform_gwas_vcf {
    tag "transform_gwas_vcf"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(gwas_vcf)
    
    output:
    path("transformed_gwas_vcf.txt.gz"), emit: transformed_gwas_vcf

    script:
    """
    echo "#CHR POS REF ALT SNP_ID BETA SE P" > temp.txt
    bcftools query -f'chr%CHROM %POS %REF %ALT [%SNP] [%BETA] [%SE] [%P]\n' $gwas_vcf >> temp.txt
    # Generating the N column
    echo "N" > n_col.txt
    for i in \$(seq 2 `wc -l < temp.txt`); do
        echo $params.gwas_sample_size >> n_col.txt
    done
    paste -d " " temp.txt n_col.txt > base.data
    calculate_z.py -i base.data -o transformed_gwas_vcf.txt
    bgzip transformed_gwas_vcf.txt
    """

}

process add_annotations {
    tag "annotate"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(vcf_sumstats)
    file(ref_fasta)
    file(ref_fasta_index) 
    file(gene_annotations) 
    file(codon_file) 
    file(priority_file) 

    output:
    path("annotated.txt.gz"), emit: annot_transformed_gwas

    script:
    """
    /anno/anno -i $vcf_sumstats -g $gene_annotations -o annotated --inputFormat plain -c $codon_file -p $priority_file -r $ref_fasta
    echo "#CHR\tPOS\tREF\tALT\tSNP_ID\tN\tZSCORE\tANNO" > annotated.txt
    tail -n +2 annotated | cut -f1-8 >> annotated.txt
    bgzip -c annotated.txt > annotated.txt.gz
    """

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

workflow lifebitai_twas_vcf{

    take:
        ch_gene_annotations
        ch_codon_file
        ch_priority_file
        ch_ref_fasta
        ch_ref_fasta_index
        ch_ld_reference
        ch_eqtl_weights
        ch_gwas_sumstats


    main:

        // Define channels from repository files

        transform_gwas_vcf(ch_gwas_sumstats)

        add_annotations(transform_gwas_vcf.out.transformed_gwas_vcf,
                    ch_ref_fasta,
                    ch_ref_fasta_index,
                    ch_gene_annotations,
                    ch_codon_file,
                    ch_priority_file)

        ptwas_scan(add_annotations.out.annot_transformed_gwas,
                    ch_ld_reference,
                    ch_eqtl_weights)

    emit:
        ptwas_scan.out.gambit_output
}

workflow {
    
    ch_gwas_sumstats = Channel
        .fromPath("${params.gwas_summary_statistics}")
        .ifEmpty { exit 1, "GWAS summary statistics file not found: ${params.gwas_summary_statistics}" }

    ch_ld_reference = Channel
    .fromPath("${params.ld_reference_panel}")
    .ifEmpty { exit 1, "File with LD reference panel not found: ${params.ld_reference_panel}" }

    ch_eqtl_weights = Channel
    .fromPath("${params.eqtl_weights}")
    .ifEmpty { exit 1, "File with eQTL weights not found: ${params.eqtl_weights}" }

    ch_gene_annotations = Channel
    .fromPath("${params.gene_annotations}")
    .ifEmpty { exit 1, "File with gene annotations not found: ${params.gene_annotations}" }

    ch_codon_file = Channel
    .fromPath("${params.codons}")
    .ifEmpty { exit 1, "File with codons not found: ${params.codons}" }

    ch_priority_file = Channel
    .fromPath("${params.priority_file}")
    .ifEmpty { exit 1, "Priority file not found: ${params.priority_file}" }

    ch_ref_fasta = Channel
        .fromPath("${params.ref_fasta}")
        .ifEmpty { exit 1, "Reference genome not found: ${params.ref_fasta}" }
    
    ch_ref_fasta_index = Channel
        .fromPath("${params.ref_fasta_index}")
        .ifEmpty { exit 1, "Reference genome index not found: ${params.ref_fasta_index}" }

    lifebitai_twas_vcf(ch_gene_annotations,
                        ch_codon_file,
                        ch_priority_file,
                        ch_ref_fasta,
                        ch_ref_fasta_index,
                        ch_ld_reference,
                        ch_eqtl_weights,
                        ch_gwas_sumstats)
}
