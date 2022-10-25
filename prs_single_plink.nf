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
target_pheno              : ${params.target_pheno}
ldpred2                   : ${params.ldpred2}
megaprs                   : ${params.megaprs}
prscsx                    : ${params.prscsx}
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

process transform_gwas_vcf_ldpred2 {
    label 'gwas_vcf'
    publishDir "${params.outdir}/transformed_VCF_input", mode: "copy"

    input:
    path gwas_vcf

    output:
    path("base.data"), emit: sumstats

    script:
    """
    echo "rsid,chr,pos,a0,a1,beta,beta_se,p,N" > base.data
    bcftools norm --multiallelic -any -Ou $gwas_vcf | \
    bcftools query -f'[%ID],%CHROM,%POS,%REF,%ALT,[%ES],[%SE],[%LP],[%SS]\n' | \
    awk 'BEGIN{FS=OFS=","} {
      print \$1,\$2,\$3,\$4,\$5,\$6,\$7,sprintf("%.10f",10^-\$8),\$9;
    }' >> base.data
    """
}

process transform_gwas_vcf_prscsx {
    label 'gwas_vcf'
    publishDir "${params.outdir}/transformed_VCF_input", mode: "copy"

    input:
    path gwas_vcf

    output:
    path("base.data"), emit: sumstats
    path("max_sample_size.txt"), emit: sample_size

    script:
    """
    printf "SNP\tA1\tA2\tBETA\tP\n" > base.data
    bcftools norm --multiallelic -any -Ou $gwas_vcf | \
    bcftools query -f'[%ID]\t%REF\t%ALT\t[%ES]\t[%LP]\n' | \
    awk 'BEGIN{FS=OFS="\t"} {
      print \$1,\$2,\$3,\$4,sprintf("%.10f", 10^-\$5);
    }' >> base.data
    bcftools query -f'[%SS]\n' $gwas_vcf | \
    awk '\$1>m{m=\$1}END{print m}' > max_sample_size.txt
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

process transform_target_pheno {

    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"

    input:
    path pheno

    output:
    tuple path("target.pheno"), path("target.cov"), emit: transformed_target_pheno //killed: transformed_target_pheno_for_plots_ch, transformed_target_pheno_for_ldpred_ch)

    script:
    """
    transform_target_pheno.R --input_pheno ${pheno}
    """

}


process estimate_per_pred_h2 {
    label 'megaprs'
    publishDir "${params.outdir}/MEGAPRS", mode: "copy"

    input:
    tuple val(plink_prefix), path(bed), path(bim), path(fam)

    output:
    path("bld65"), emit: bld65

    script:
    """
    ldak5.2.linux --cut-weights sections --bfile $plink_prefix
    ldak5.2.linux --calc-weights-all sections --bfile $plink_prefix
    mv sections/weights.short bld65
    """

}

process calc_per_pred_hers {
    label 'megaprs'
    publishDir "${params.outdir}/MEGAPRS", mode: "copy"

    input:
    tuple val(plink_prefix), path(bed), path(bim), path(fam)
    path(pheno)
    path(bld65)
    path('*')

    output:
    path("bld.ldak.ind.hers"), emit: per_pred_hers
    path("quant.summaries"), emit: sumstats_megaprs

    script:
    """
    ldak5.2.linux --linear quant --bfile $plink_prefix --pheno $pheno
    ldak5.2.linux --calc-tagging bld.ldak \
                            --bfile $plink_prefix \
                            --ignore-weights YES \
                            --power -.25 \
                            --annotation-number 65 \
                            --annotation-prefix bld \
                            --window-cm 1 \
                            --save-matrix YES
    ldak5.2.linux --sum-hers bld.ldak \
                    --tagfile bld.ldak.tagging \
                    --summary quant.summaries \
                    --matrix bld.ldak.matrix
    """

}

process detect_high_ld_snps { 

    label 'megaprs'
    publishDir "${params.outdir}/MEGAPRS", mode: "copy"

    input:
    tuple val(plink_prefix), path(bed), path(bim), path(fam)
    path(pheno)
    path(highld)
    tuple val(ref_plink_prefix), path(ref_bed), path(ref_bim), path(ref_fam)

    output:
    path("cors*"), emit: pred_cors
    path("highld/genes.predictors.used"), emit: high_ld_pred


    script: 
    """
    ldak5.2.linux --cut-genes highld --bfile $ref_plink_prefix --genefile $highld
    ldak5.2.linux --calc-cors cors --bfile $plink_prefix --window-cm 3
    chmod 777 highld/
    """

    }

process megaprs_pred_model {

    label 'megaprs'
    publishDir "${params.outdir}/MEGAPRS", mode: "copy"

    input:
    path("*")
    path(ss)
    path(ldak_ind_hers)
    path(high_ld_pred)

    output:
    path("megaprs.effects"), emit: megaprs_effect_weights

    script:
    """
    ldak5.2.linux --mega-prs megaprs --model mega \
                        --ind-hers $ldak_ind_hers \
                        --summary $ss \
                        --cors cors \
                        --cv-proportion .1 \
                        --high-LD $high_ld_pred \
                        --window-cm 1 \
                        --allow-ambiguous YES
    """
    
}

process prepare_weights_file {

    publishDir "${params.outdir}/MEGAPRS", mode: "copy"
        
    input:
    path(weights_file)
    tuple val(plink_prefix), path(bed), path(bim), path(fam)
        
    output:
    path("clean_megaprs*"), emit: megaprs_filled_weights

    script:
    """
    prepare_weights.R \
                --variant_weights=${weights_file} \
                --bim_file=${bim} \
                --output=clean_megaprs

    """
}

process calculate_sample_prs {
    label "plink2"
    publishDir "${params.outdir}/MEGAPRS/PRS", mode: "copy"
        
    input:  
    tuple val(plink_prefix), path(bed), path(bim), path(fam)   
    path(clean_megaprs_weights)

    output:
    path("megaprs_prs.sscore"), emit: megaprs_prs
    path("megaprs_prs_harmonised.sscore"), emit: scores_to_percentiles

    script:
    """
    plink2 --bfile $plink_prefix \
        --score $clean_megaprs_weights header-read \
        --score-col-nums 3 \
        --out megaprs_prs
    # Create a suitable version of the scores file for PRS percentile calculation
    cp megaprs_prs.sscore megaprs_prs_harmonised.sscore
    sed -i -E '1 s/SCORE_.+_AVG/PRS/' megaprs_prs_harmonised.sscore
    """
}



process run_ldpred2 {
    label "ldpred2"
    publishDir "${params.outdir}/LDpred2", mode: "copy"

    input:
    tuple val(plink_genotypes_name), path(bed), path(bim), path(fam)
    path(summary_statistics)
    path(ld_reference_rds)

    output: 
    path("ldpred2_prs_scores.csv"), emit: ldpred2_prs_scores
    path("ldpred2_scores_harmonised.csv"), emit: scores_to_percentiles
    
    script:
    """
    ldpred2.R --plink_prefix=$plink_genotypes_name \
                                --gt_format=${params.gt_format} \
                                --summary_statistics=$summary_statistics \
                                --validation_test_split=${params.validation_test_split} \
                                --ld_reference=$ld_reference_rds \
                                --random_seed_split=${params.random_seed_split} \
                                --output_file=ldpred2_prs_scores.csv
    # Create a suitable version of the scores file for PRS percentile calculation
    cp ldpred2_prs_scores.csv ldpred2_scores_harmonised.csv
    sed -i 's/beta_auto/PRS/' ldpred2_scores_harmonised.csv
    """
}

process estimate_posterior_effects_prscsx {
    label 'prscsx'
    publishDir "${params.outdir}/PRS-CSX", mode: "copy"

    input:
    tuple val(plink_prefix), path(bed), path(bim), path(fam)
    path(sumstats_file)
    path(ldref)
    path(snpinfo)
    path(sample_size)

    output:
    path("results/*"), emit: output_prscsx

    script:
    """
    #create ref dir
    tar -xzvf ${ldref}
    mkdir reference
    mkdir data
    mkdir results
    mv ${ldref.simpleName} reference/
    cp ${snpinfo} reference/
    cp ${bim} data/
    cp ${sumstats_file} data/
    chr=${params.chromosome}
    gwas_sample_size=\$(grep -E '[0-9]' ${sample_size})
    if [ \$chr = false ]; then
        python /PRScsx/PRScsx.py \
            --ref_dir=reference \
            --bim_prefix=data/${bim.simpleName} \
            --sst_file=data/${sumstats_file} \
            --n_gwas=\$gwas_sample_size \
            --pop=${params.population} \
            --phi=${params.phi} \
            --out_dir=results \
            --out_name=${params.outfile_name}
    else
        python /PRScsx/PRScsx.py \
            --ref_dir=reference \
            --bim_prefix=data/${bim.simpleName} \
            --sst_file=data/${sumstats_file} \
            --n_gwas=\$gwas_sample_size \
            --pop=${params.population} \
            --chrom=${params.chromosome} \
            --phi=${params.phi} \
            --out_dir=results \
            --out_name=${params.outfile_name}
    fi
    rm *.tar.gz
    rm -r reference
    """

}

process calculate_sample_prs_prscsx {
    label "plink2"
    publishDir "${params.outdir}/PRS-CSX/results/PRS", mode: "copy"
        
    input:  
    tuple val(plink_prefix), path(bed), path(bim), path(fam)     
    path("results/*")

    output:
    path("prscsx_prs.sscore")
    path("prscsx_prs_harmonised.sscore"), emit: scores_to_percentiles

    script:
    """
    nfiles=`ls results/*.txt|wc -l`
    if [ \$nfiles -gt 1 ]; then
        cat results/*.txt > results/scores.txt
    fi
    if [ \$nfiles -eq 1 ]; then
        mv results/*.txt results/scores.txt
    fi
    plink2 --bfile $plink_prefix \
        --score results/scores.txt 2 4 6 cols=+scoresums \
        --out prscsx_prs
    # Create a suitable version of the scores file for PRS percentile calculation
    cp prscsx_prs.sscore prscsx_prs_harmonised.sscore
    sed -i 's/SCORE1_SUM/PRS/' prscsx_prs_harmonised.sscore
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

process ld_prune {
    label "plink"

    input:
    tuple val(plink_genotypes_name), path(bed), path(bim), path(fam)

    output:
    tuple val('ldpruned'), path('ldpruned.bed'), path('ldpruned.bim'), path('ldpruned.fam'), emit: pruned_snps //killed: ch_pruned_snps_for_pca, ch_pruned_snps_for_grm)

    script:
    """
    plink \
        --bfile ${plink_genotypes_name} \
        --keep-allele-order \
        --allow-no-sex \
        --geno ${params.grm_miss} \
        --maf ${params.grm_maf} \
        --mac ${params.grm_mac} \
        --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
        --memory ${task.memory.toMega()} \
        --out filtered_ldpruned

    plink \
        --bfile ${plink_genotypes_name} \
        --keep-allele-order \
        --allow-no-sex \
        --extract filtered_ldpruned.prune.in \
        --make-bed \
        --memory ${task.memory.toMega()} \
        --out ldpruned
    """
}

process calculate_pcs {
    label "plink"

    input:
    tuple val(name), path('ldpruned.bed'), path('ldpruned.bim'), path('ldpruned.fam')

    output:
    tuple path('pca.eigenvec'), path('pca.eigenval'), emit: pc_table

    script:
    """
    plink \
        --bfile ldpruned \
        --keep-allele-order \
        --pca ${params.number_pcs} header tabs \
        --memory ${task.memory.toMega()} \
        --out pca
    """
}

process calculate_grm {
    label "plink"

    input:
    tuple val(name), path('ldpruned.bed'), path('ldpruned.bim'), path('ldpruned.fam')

    output:
    tuple path('grm.grm.bin'), path('grm.grm.N.bin'), path('grm.grm.id'), emit: grm_out

    script:
    """
    plink \
        --bfile ldpruned \
        --keep-allele-order \
        --make-grm-bin \
        --memory ${task.memory.toMega()} \
        --out grm
    """
}

process create_post_analysis_dir {
    label "prs_notebook"
    publishDir "${params.outdir}/combined_results", mode: "copy", 
        saveAs: {filename -> filename.split('/')[-1]}

    input:
    tuple path('grm.grm.bin'), path('grm.grm.N.bin'), path('grm.grm.id')
    tuple path('pca.eigenvec'), path('pca.eigenval')
    path("phenofile.phe")
    path("*")
    path('PRS_post_analysis_norun.ipynb')

    output:
    path("combined/*")

    script:
    """
    mkdir combined && find -L . -type f -not -path '*/.*'  -not -path '*/*.ipynb' | xargs -I {} cp {} combined
    sed 's/<phenocolplaceholder>/${params.pheno_col}/' PRS_post_analysis_norun.ipynb > combined/PRS_post_analysis.ipynb
    cd combined
    jupyter nbconvert --to notebook --execute PRS_post_analysis.ipynb \
        --allow-errors \
        --inplace \
        --ExecutePreprocessor.kernel_name="ir"
    """
}


workflow lifebitai_prs_single_plink{
        take:
            ch_gwas_vcf
   
        main:
            ch_prs_scores_tables = Channel.empty()

            ch_genotypes = Channel
                .fromFilePairs("${params.genotype_data}",size:3, flat : true)
                .ifEmpty { exit 1, "Genotype data in plink format not found: ${params.genotype_data}" }
            
            if (params.target_pheno) {
                ch_target_pheno = Channel
                    .fromPath(params.target_pheno, checkIfExists: true)
                    .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
            }

            if (params.ldpred2) {
                if (params.ld_reference_rds) {
                    ch_ld_reference_rds = Channel
                        .fromPath(params.ld_reference_rds, checkIfExists: true)
                        .ifEmpty { exit 1, "LD reference .rds file not found: ${params.ld_reference_rds}" }
                }
            }
            if (params.megaprs) {
                if (params.target_pheno) {
                    ch_pheno_file = Channel
                        .fromPath(params.target_pheno, checkIfExists: true)
                        .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
                }
                if (params.high_ld_file) {
                    ch_highld_file = Channel
                                            .fromPath("${params.high_ld_file}")
                                            .ifEmpty { exit 1, "File with high LD regions not found: ${params.high_ld_file}" }
                }
                if (params.bld_annotations) {
                    ch_bld_annotations = Channel
                                                .fromPath("${params.bld_annotations}")
                                                .ifEmpty { exit 1, "BLD-LDAK annotations not found: ${params.bld_annotations}" }
                                                .collect()
                }
            }

            if (params.megaprs) {
                if (params.reference_panel) {
                    ch_ref = Channel
                        .fromFilePairs("${params.reference_panel}",size:3, flat : true)
                        .ifEmpty { exit 1, "Reference files in plink format not found: ${params.reference_panel}" }
                }
            }

            if (params.prscsx) {
                ch_ld_ref_file = !params.ld_ref_file && params.population ? Channel.value(file(params.pops[params.population].ld_ref_file)) :  Channel.value(file(params.ld_ref_file))
                ch_snpinfo_file = !params.snpinfo_file && params.population ? Channel.value(file(params.pops[params.population].snpinfo_file)) :  Channel.value(file(params.snpinfo_file))
            }
            /*----------------------------
            Setting up other parameters
            ------------------------------*/

            // // Initialise variable to store optional parameters
            if (params.ldpred2) {
                transform_gwas_vcf_ldpred2(ch_gwas_vcf)
                ch_input_sumstats = transform_gwas_vcf_ldpred2.out.sumstats
            }

            if (params.prscsx) {
                transform_gwas_vcf_prscsx(ch_gwas_vcf)
                ch_input_sumstats = transform_gwas_vcf_prscsx.out.sumstats
            }


            if ( params.prscsx) {
                // PRS-CSX improves cross-population polygenic prediction by integrating GWAS summary statistics from multiple populations. 
                // Step 1: Calculates posterior SNP effect size estimates for each chromosome.
                // Step 2: Calculates individual-level polygenic score.
                estimate_posterior_effects_prscsx(ch_genotypes,
                                                    ch_input_sumstats,
                                                    ch_ld_ref_file,
                                                    ch_snpinfo_file,
                                                    transform_gwas_vcf_prscsx.out.sample_size)
                
                calculate_sample_prs_prscsx(ch_genotypes,
                                            estimate_posterior_effects_prscsx.out.output_prscsx)
                
                ch_scores_to_percentiles = calculate_sample_prs_prscsx.out.scores_to_percentiles
            }

            if ( params.megaprs) {
                // Estimate per-predictor heritabilities assuming the BLD-LDAK Model (recommended for human data)

                // Step 1: Cut predictors into sections
                // Step 2: Calculates LDAK weightings and joins them up

                estimate_per_pred_h2(ch_genotypes)

                calc_per_pred_hers(ch_genotypes,
                                    ch_pheno_file,
                                    estimate_per_pred_h2.out.bld65,
                                    ch_bld_annotations)
                
                detect_high_ld_snps(ch_genotypes,
                                    ch_pheno_file,
                                    ch_highld_file,
                                    ch_ref)
                
                megaprs_pred_model(detect_high_ld_snps.out.pred_cors,
                                    calc_per_pred_hers.out.sumstats_megaprs,
                                    calc_per_pred_hers.out.per_pred_hers,
                                    detect_high_ld_snps.out.high_ld_pred)

                prepare_weights_file(megaprs_pred_model.out.megaprs_effect_weights,
                                        ch_genotypes)
                
                calculate_sample_prs(ch_genotypes,
                                        prepare_weights_file.out.megaprs_filled_weights)
                
                ch_prs_scores_tables = ch_prs_scores_tables.mix(calculate_sample_prs.out.megaprs_prs)

                ch_scores_to_percentiles = calculate_sample_prs.out.scores_to_percentiles

            }

            if ( params.ldpred2) {
                run_ldpred2(ch_genotypes,
                            ch_input_sumstats,
                            ch_ld_reference_rds)
                
                ch_scores_to_percentiles = run_ldpred2.out.scores_to_percentiles
                ch_prs_scores_tables = ch_prs_scores_tables.mix(run_ldpred2.out.ldpred2_prs_scores)
            }

            if (params.calc_prs_percentiles) {
            calculate_prs_percentiles(ch_scores_to_percentiles)
            }

            if (params.filter_by_percentiles) {

                filter_by_percentiles_v2(ch_genotypes,
                                        calculate_prs_percentiles.out.samples_filtered_by_percentile)
                
                ch_filtered_percentile_cohort = filter_by_percentiles_v2.out.filtered_percentile_cohort
            }

            if (params.post_analysis) {
                ld_prune(ch_genotypes)

                calculate_pcs(ld_prune.out.pruned_snps)

                calculate_grm(ld_prune.out.pruned_snps)

                post_analysis_notebook = Channel.fromPath("$projectDir/bin/PRS_post_analysis.ipynb")


                create_post_analysis_dir(calculate_grm.out.grm_out,
                                        calculate_pcs.out.pc_table,
                                        ch_target_pheno,
                                        ch_prs_scores_tables.collect(),
                                        post_analysis_notebook)

            }
        
        emit:
            ch_prs_scores_tables
}

workflow{
    allowed_gt_list = allowed_genotypic_formats()
    // Check inputs

    // Check that one and only one tool is called at a time
    tools_called = (params.ldpred2 ? 1 : 0) + \
                (params.megaprs ? 1 : 0) + \
                (params.prscsx ? 1 : 0)
    if (tools_called != 1) {
        exit 1, """Please select one, and only one, of the included tools to run:
                --ldpred2
                --megaprs
                --prscsx"""
    }

    if (params.post_analysis && !(params.ldpred2 || params.megaprs)) {
        exit 1, "Post analysis only supports LDpred2 and MEGAPRS. Please use either --ldpred2 or --megaprs."
    }

    if (params.post_analysis && !params.target_pheno) {
        exit 1, "Post analysis requires a phenofile. Please supply one using --target_pheno."
    }

    /*-----------------------------------------------------------------------------------------
    Setting up target dataset: PLINK files and pheno file output from lifebit-ai/biobank-gwas
    --------------------------------------------------------------------------------------------*/

    projectDir = workflow.projectDir

    ch_gwas_vcf = Channel
                    .fromPath(params.gwas_vcf, checkIfExists: true)
                    .ifEmpty { exit 1, "GWAS summary stats VCF not found: ${params.gwas_vcf}" }

    lifebitai_prs_single_plink(ch_gwas_vcf,
                                ch_genotypes)

}
