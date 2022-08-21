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
    python3 calculate_z.py -i base.data -o transformed_gwas_vcf.txt
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

workflow lifebitaixf_twas_vcf{

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
        ch_ptwas_gwas_sumstats = params.annotate_vcf ? transform_gwas_vcf.out.transformed_gwas_vcf : ch_transformed_gwas_vcf

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
