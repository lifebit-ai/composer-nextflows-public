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
genotype_data             : ${params.genotype_data}
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

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Define list of available tools to annotate
def allowed_genotypic_formats() {
    return [
        'dosages',
        'hard_called'
    ]
}


process retrieve_pgs_weights {
    label "python"
    publishDir "${params.outdir}/PGS_catalog", mode: "copy"

    input:
    val(trait_id)

    output:
    path("*.txt.gz"), emit: pgs_weights
    path("*.json"), emit: pgs_metadata_json

    script:
    """
    curl -X GET 'https://www.pgscatalog.org/rest/score/search?trait_id=${trait_id}' -H  "accept: application/json" > ${trait_id}_query_result.json
    process_json.py -j ${trait_id}_query_result.json
    """
}

process prepare_pgs_weights_file {
    publishDir "${params.outdir}/PGS_catalog", mode: "copy"

    input:
    path(traits)

    output:
    path("pgs.score"), emit: pgs_weights_clean

    script:
    """
    echo "SNP A1 SCORE_PGS" > pgs.score
    # Check the first column
    FIRST_COL=\$( zcat $traits | sed '/^#/d' | head -1 | awk -F "\t" '{print \$1}')
    if [ \$FIRST_COL == "rsID" ]; then
        zcat $traits | sed '/^#/d' | sed '1d' | awk 'BEGIN{OFS=" ";FS="\t"}{print \$2":"\$3,\$4,\$6}' >> pgs.score
    else
        zcat $traits | sed '/^#/d' | sed '1d' | awk 'BEGIN{OFS=" ";FS="\t"}{print \$1":"\$2,\$3,\$5}' >> pgs.score
    fi
    """
}

process calculate_sample_prs_with_pgs_weights {
    label "plink2"
    publishDir "${params.outdir}/PGS_catalog/PRS", mode: "copy"

    input:  
    tuple val(plink_prefix), path(bed), path(bim), path(fam)    
    path(clean_pgs_weights)

    output:
    path("PGS_prs.sscore")
    path("PGS_prs_harmonised.sscore"), emit: scores_to_percentiles

    script:
    """
    plink2 --bfile $plink_prefix \
        --score $clean_pgs_weights header-read \
        --score-col-nums 3 \
        --out PGS_prs
    # Create a suitable version of the scores file for PRS percentile calculation
    cp PGS_prs.sscore PGS_prs_harmonised.sscore
    sed -i -E '1 s/SCORE_.+_AVG/PRS/' PGS_prs_harmonised.sscore
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

workflow lifebitai_prs_pgs{
        take:
            ch_pgs_catalogue_weights

        main:
            ch_genotypes = Channel
                                .fromFilePairs("${params.genotype_data}",size:3, flat : true)
                                .ifEmpty { exit 1, "Genotype data in plink format not found: ${params.genotype_data}" }

            retrieve_pgs_weights(ch_pgs_catalogue_weights)
    
            prepare_pgs_weights_file(retrieve_pgs_weights.out.pgs_weights)

            calculate_sample_prs_with_pgs_weights(ch_genotypes,
                                                    prepare_pgs_weights_file.out.pgs_weights_clean)
            
            ch_scores_to_percentiles = calculate_sample_prs_with_pgs_weights.out.scores_to_percentiles

            if (params.calc_prs_percentiles) {
                calculate_prs_percentiles(ch_scores_to_percentiles)
            }

            if (params.filter_by_percentiles) {
                filter_by_percentiles_v2(ch_genotypes,
                                        calculate_prs_percentiles.out.samples_filtered_by_percentile)
                
                ch_filtered_percentile_cohort = filter_by_percentiles_v2.out.filtered_percentile_cohort
            }
            
        
        emit:
            ch_scores_to_percentiles
            
}

workflow{

    ch_pgs_catalogue_weights = Channel
        .from(params.pgs_weights_by_trait_id)
        .ifEmpty { exit 1, "trait_id to retrieve PGS catalogue weight not found: ${params.pgs_weights_by_trait_id}"}

    lifebitai_prs_pgs(ch_pgs_catalogue_weights)

}