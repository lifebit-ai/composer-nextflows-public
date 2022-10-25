#!/usr/bin/env nextflow
/*
===============================================================================
                         bi-polygenic-risk-scores-nf
===============================================================================
 Polygenic Risk Scores Pipeline to calculate polygenic risk scores.

-------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

log.info """\
P R S  P I P E L I N E
==========================================================
gwas_vcf                  : ${params.gwas_vcf}
calc_prs_percentiles      : ${params.calc_prs_percentiles}
filter_by_percentiles     : ${params.filter_by_percentiles}
outdir                    : ${params.outdir}             
"""

/*-----------------------------------------------------------------------------------------
  Functions for validating parameters and inputs
--------------------------------------------------------------------------------------------*/

def get_chromosome( file ) {
    // using RegEx to extract chromosome number from file name
    regexpPE = /(?:chr)[a-zA-Z0-9]+/
    (file =~ regexpPE)[0].replaceAll('chr','').toInteger()
}


process transform_gwas_vcf_sbayesr {
    label 'gwas_vcf'
    publishDir "${params.outdir}/transformed_VCF_input", mode: "copy"

    input:
    path gwas_vcf

    output:
    path("base.data"), emit: sumstats

    script:
    """
    echo "SNP,CHR,POS,A1,A2,FREQ,BETA,SE,P,N" > base.data
    bcftools norm --multiallelic -any -Ou $gwas_vcf | \
    bcftools query -f'[%ID],%CHROM,%POS,%REF,%ALT,[%AF],[%ES],[%SE],[%LP],[%SS]\n' | \
    awk 'BEGIN{FS=OFS=","} {
      print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,sprintf("%.10f",10^-\$9),\$10;
    }' >> base.data
    """
}

process split_sumstats_sbayesr {
    label "high_memory"
    publishDir "${params.outdir}/sbayesr", mode: "copy"

    input:
    each path(sumstats)
    val(chrom)

    output:
    tuple val(chrom), path("GWAS_sumstats_COJO.chr*.txt"), emit: split_sumstats_output

    script:
    """
    split_sumstats_sbayesr.R --input_sumstats=$sumstats \
                                            --chromosome=$chrom \
                                            --p_max=${params.p_max}
    """
}

process run_sbayesr {
    label "sbayesr"
    publishDir "${params.outdir}/sbayesr", mode: "copy"
    stageInMode 'copy'
    
    input:
    tuple val(chrom), path(sumstats), path(bin), path(info)

    output:
    path("*.snpRes"), emit: sbayesr_snp_weights
    script:
    """
    gctb  --sbayes R \
    --ldm ${bin.baseName} \
    --pi ${params.pi_values} \
    --gamma ${params.gamma_values} \
    --gwas-summary ${sumstats} \
    --robust --exclude-mhc --out-freq 10 \
    --out sbayesr_chr${chrom} \
    --chain-length 10000 --burn-in 2000 2>&1 | tee 'sim_${chrom}.log'
    """
}


process concat_sbayesr_weights {
    label "sbayesr"
    publishDir "${params.outdir}/sbayesr", mode: "copy"
    stageInMode 'copy'
    
    input:
    path(chr_weights)

    output:
    path("sbayesR_gw_weights.score.gz"), emit: gw_sbayesr_weights

    script:
    """
    awk 'FNR==1 && NR!=1 { while (/^Id/) getline; }1 {print}' *.snpRes > gw_weights.txt

    echo "SNP A1 SCORE_SBayesR" > sbayesR_gw_weights.score
    awk '{print \$2,\$5,\$8}' gw_weights.txt >> sbayesR_gw_weights.score
    gzip sbayesR_gw_weights.score

    """
}

process calculate_sample_prs_sbayesr {
    label "plink2"
    publishDir "${params.outdir}/sbayesr/PRS", mode: "copy"
        
    input:  
    tuple val(plink_prefix), path(bed), path(bim), path(fam)   
    path(gw_sbayesr_weights)

    output:
    path("sbayesr_prs.sscore"), emit: sbayesr_prs
    path("sbayesr_prs_harmonised.sscore"), emit: scores_to_percentiles

    script:
    """
    plink2 --bfile $plink_prefix \
        --score $gw_sbayesr_weights header-read \
        --score-col-nums 3 \
        --out sbayesr_prs
    # Create a suitable version of the scores file for PRS percentile calculation
    cp sbayesr_prs.sscore sbayesr_prs_harmonised.sscore
    sed -i -E '1 s/SCORE_.+_AVG/PRS/' sbayesr_prs_harmonised.sscore
    """
}

process calculate_prs_percentiles {
    label "python"
    publishDir "${params.outdir}/prs_percentiles", mode: "copy"

    input:
    path(prs_scores)

    output:
    path("prs_scores_percentiles.csv"), emit: prs_scores_percentiles
    path("samples_filtered_percentile.tsv"), emit: samples_filtered_by_percentile

    script:
    """
    calculate_prs_percentiles.py --input_prs_scores $prs_scores \
                                    --output_prs_df prs_scores_percentiles.csv \
                                    --prs_column_name ${params.prs_column_name} \
                                    --lower_prs_percentile ${params.lower_prs_percentile} \
                                    --upper_prs_percentile ${params.upper_prs_percentile} 

    """
}

process merge_plink {

    label "plink"
    publishDir "${params.outdir}/merged_plink", mode: "copy"

    input:
    tuple val(name), path("*")

    output:
    path("merged.*"), emit: merged_plink

    script:
    """ 
    ls *.bed > bed.txt
    ls *.bim > bim.txt
    ls *.fam > fam.txt

    FIRSTFILE=\$(head -1 bed.txt > first_file.txt && while IFS= read -r line; do echo \${line%%.*}; done < first_file.txt)

    sed -i '1d' bed.txt
    sed -i '1d' bim.txt
    sed -i '1d' fam.txt

    paste bed.txt bim.txt fam.txt > merge.list

    plink --keep-allele-order --bfile \$FIRSTFILE --merge-list merge.list --allow-no-sex --make-bed --out merged
    """
}

process filter_by_percentiles_v1 {
    label "plink"
    publishDir "${params.outdir}/filtered_${params.lower_prs_percentile}_${params.upper_prs_percentile}_percentile", mode: "copy"

    input:
    each path(merged_plink_file)
    path(samples_filtered_by_percentile)

    output:
    path("filtered_by_*"), emit: filtered_percentile_cohort

    script:
    """
    plink --bfile merged \
    --keep-allele-order \
    --keep $samples_filtered_by_percentile \
    --make-bed \
    --out filtered_by_${params.lower_prs_percentile}_${params.upper_prs_percentile}_percentile
    """
}


process filter_by_percentiles_v2 {
    label "plink"
    publishDir "${params.outdir}/filtered_${params.lower_prs_percentile}_${params.upper_prs_percentile}_percentile", mode: "copy"

    input:
    tuple val(plink_prefix), path(bed), path(bim), path(fam)
    path(samples_filtered_by_percentile)

    output:
    path("filtered_by_*"), emit: filtered_percentile_cohort

    script:
    """
    # Reformat fam file to have FID and IID with the same values as required
    cp $fam temp.fam
    cut --complement -d "\t" -f1 $fam > temp.fam
    awk -F "\t" '{print \$1}' temp.fam > col1.txt
    paste -d "\t" col1.txt temp.fam > ${plink_prefix}.fam
    rm temp.fam col1.txt
    plink --bfile $plink_prefix \
    --keep-allele-order \
    --keep $samples_filtered_by_percentile \
    --make-bed \
    --out filtered_by_${params.lower_prs_percentile}_${params.upper_prs_percentile}_percentile
    """
}

workflow lifebitai_prs_sbayesr{
    take:
        ch_gwas_vcf

    main:
        transform_gwas_vcf_sbayesr(ch_gwas_vcf)

        ch_ref = Channel
            .fromFilePairs("${params.reference_panel}",size:3, flat : true)
            .ifEmpty { exit 1, "Reference files in plink format not found: ${params.reference_panel}" }

        if ( params.chrom_range) {
            ch_chromosomes = Channel
                .from( params.chrom_range )
                .ifEmpty { exit 1, "Chromosome range not supplied: ${params.chrom_range}" }
        }

        if ( params.chr_ld_matrix_dir ) {
            ch_chr_ld_matrix = Channel
                .fromPath("${params.chr_ld_matrix_dir}/*chr${params.chrom_range}*.{bin,info}" )
                .map { it -> [ get_chromosome(file(it).simpleName.minus(".ldm.sparse.bin").minus(".ldm.sparse.info")), "s3:/"+it] }
                .groupTuple(by:0,sort: true)
                .map { chr, files_pair -> [ chr, files_pair[0], files_pair[1] ] }
                .map { chr, bin, info -> [ chr, file(bin), file(info) ] }
        }

        split_sumstats_sbayesr(transform_gwas_vcf_sbayesr.out.sumstats,
                                ch_chromosomes)
        
        ch_sbayesr_input = split_sumstats_sbayesr.out.split_sumstats_output
                            .join(ch_chr_ld_matrix,by: 0) // killed: ch_sbayesr_test
        
        run_sbayesr(ch_sbayesr_input)

        concat_sbayesr_weights(run_sbayesr.out.sbayesr_snp_weights.collect())


        calculate_sample_prs_sbayesr(ch_ref,
                                        concat_sbayesr_weights.out.gw_sbayesr_weights)
        
        ch_scores_to_percentiles = calculate_sample_prs_sbayesr.out.scores_to_percentiles

        if (params.calc_prs_percentiles) {
            calculate_prs_percentiles(ch_scores_to_percentiles)
        }

        if (params.filter_by_percentiles) {
            if (params.merge_per_chrom) {
                filter_by_percentiles_v1(merge_plink.out.merged_plink,
                                            calculate_prs_percentiles.out.samples_filtered_by_percentile)
            
                ch_filtered_percentile_cohort = filter_by_percentiles_v1.out.filtered_percentile_cohort
            
            } else {
                filter_by_percentiles_v2(ch_ref,
                                        calculate_prs_percentiles.out.samples_filtered_by_percentile)
                
                ch_filtered_percentile_cohort = filter_by_percentiles_v2.out.filtered_percentile_cohort
            }
        }
    
    emit:
        ch_scores_to_percentiles
}

workflow{
    // Check inputs

    ch_prs_scores_tables = Channel.empty()

    projectDir = workflow.projectDir


    if (params.gwas_vcf) {
        ch_gwas_vcf = Channel
        .fromPath(params.gwas_vcf, checkIfExists: true)
        .ifEmpty { exit 1, "GWAS summary stats VCF not found: ${params.gwas_vcf}" }
    }

    /*----------------------------
    Setting up other parameters
    ------------------------------*/

    if ( params.sbayesr) {
       lifebitai_prs_sbayesr(
            ch_gwas_vcf
       )
    }
}
