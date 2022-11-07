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
ldpred1                   : ${params.ldpred1}
prsice                    : ${params.prsice}
genotype_data             : ${params.genotype_data}
pheno_metadata            : ${params.pheno_metadata}
target_plink_files_dir    : ${params.target_plink_dir}
target_pheno              : ${params.target_pheno}
binary_trait              : ${params.binary_trait}
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

/*---------------------------------------
  Transforming GWAS VCF input MRCIEU v1.0
----------------------------------------*/
process transform_saige_base {
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"

    input:
    path saige_base

    output:
    path("base.data"), emit: transformed_base //killed: transformed_base_ldpred_ch)
        
    script:
    """
    transform_base_saige.R --input_saige ${saige_base}
    """
}

process transform_gwas_vcf_ldpred1 {
    label 'gwas_vcf'
    publishDir "${params.outdir}/transformed_VCF_input", mode: "copy"

    input:
    path gwas_vcf

    output:
    path("base.data"), emit: transformed_base_ldpred

    script:
    """
    echo "CHR POS SNPID Allele1 Allele2 BETA SE p.value N" > base.data
    bcftools norm --multiallelic -any -Ou $gwas_vcf | \
    bcftools query -f'%CHROM %POS [%ID] %REF %ALT [%ES] [%SE] [%LP] [%SS]\n' | \
    awk 'BEGIN{FS=OFS=" "} {
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,sprintf("%.10f", 10^-\$8),\$9;
    }' >> base.data
    """
}

process transform_gwas_vcf_prsice {
    label 'gwas_vcf'
    publishDir "${params.outdir}/transformed_VCF_input", mode: "copy"

    input:
    path gwas_vcf

    output:
    path("base.data"), emit: transformed_base

    script:
    """
    echo "CHR POS SNPID Allele1 Allele2 BETA SE p.value" > base.data.pre
    bcftools norm --multiallelic -any -Ou ${gwas_vcf} | \
    bcftools query -f'%CHROM %POS [%ID] %REF %ALT [%ES] [%SE] [%LP]\n' | \
    awk 'BEGIN{FS=OFS=" "} {
      print \$1,\$2,\$3,\$4,\$5,\$6,\$7,sprintf("%.10f", 10^-\$8);
    }' >> base.data.pre
    # Remove duplicate SNPIDs to deal with PRSice unable to handle multiallelic sites
    awk '!seen[\$3]++' base.data.pre > base.data
    """
}

process merge_plink {

    label "plink"
    publishDir "${params.outdir}/merged_plink", mode: "copy"

    input:
    tuple val(name), path("*")

    output:
    path("merged.*"), emit: merged_plink //killed: merged_plink_ldpred2_ch, merged_plink_ldpred_gibbs_ch, merged_plink_ldpred_score_ch, merged_plink_percentile )

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

process polygen_risk_calcs {

    publishDir "${params.outdir}", mode: "copy"

    input:
    path base
    tuple val(name), path("*")
    tuple path(pheno), path(cov)

    output:
    path("*"), emit: all_results
    path("PRSice.best"), emit: best_PRS //killed: scores_to_percentiles_ch)
    path("PRSice*"), emit: results_for_report

    shell:
    extra_flags = ""

    // 1 - Base file options

    if ( params.base_info ) { extra_flags += " --base-info ${params.base_info}" }
    if ( params.base_maf ) { extra_flags += " --base-maf ${params.base_maf}" }
    if ( params.no_default ) { extra_flags += " --no-default ${params.no_default}" }

    // 2 - Target file options (other than the hardcoded ones in main.nf)

    if ( params.geno ) { extra_flags += " --geno ${params.geno}" }
    if ( params.info ) { extra_flags += " --info ${params.info}" }
    if ( params.keep ) { extra_flags += " --keep ${params.keep}" }
    if ( params.maf ) { extra_flags += " --maf ${params.maf}" }
    if ( params.nonfounders ) { extra_flags += " --nonfounders ${params.nonfounders}" }
    if ( params.ignore_fid ) { extra_flags += " --ignore-fid ${params.ignore_fid}" }
    if ( params.prevalence ) { extra_flags += " --prevalence ${params.prevalence}" }
    if ( params.remove ) { extra_flags += " --remove ${params.remove}" }

    // 3 - Dosage commands not available (pipeline is not currently developed to handle dosages)

    // 4 - Clumping

    if ( params.clump_kb ) { extra_flags += " --clump-kb ${params.clump_kb}" }
    if ( params.clump_r2 ) { extra_flags += " --clump-r2 ${params.clump_r2}" }
    if ( params.clump_p ) { extra_flags += " --clump-p ${params.clump_p}" }
    if ( params.ld ) { extra_flags += " --ld ${params.ld}" }
    if ( params.ld_dose_thres ) { extra_flags += " --ld-dose-thres ${params.ld_dose_thres}" }
    if ( params.ld_geno ) { extra_flags += " --ld-geno ${params.ld_geno}" }
    if ( params.ld_info ) { extra_flags += " --ld-info ${params.ld_info}" }
    if ( params.ld_hard_thres ) { extra_flags += " --ld-hard-thres ${params.ld_hard_thres}" }
    if ( params.ld_keep ) { extra_flags += " --ld-keep ${params.ld_keep}" }
    if ( params.ld_list ) { extra_flags += " --ld-list ${params.ld_list}" }
    if ( params.ld_maf ) { extra_flags += " --ld-maf ${params.ld_maf}" }
    if ( params.ld_remove ) { extra_flags += " --ld-remove ${params.ld_remove}" }
    if ( params.ld_type ) { extra_flags += " --ld-type ${params.ld_type}" }
    if ( params.no_clump ) { extra_flags += " --no-clump ${params.no_clump}" }
    if ( params.proxy ) { extra_flags += " --proxy ${params.proxy}" }

    // 5 - Covariate commands not available as pipeline handles the covariates directly

    // 6 - P-value thresholding

    if ( params.bar_levels ) { extra_flags += " --bar-levels ${params.bar_levels}" }
    if ( params.fastscore ) { extra_flags += " --fastscore ${params.fastscore}" }
    if ( params.no_full ) { extra_flags += " --no-full ${params.no_full}" }
    if ( params.interval ) { extra_flags += " --interval ${params.interval}" }
    if ( params.lower ) { extra_flags += " --lower ${params.lower}" }
    if ( params.model ) { extra_flags += " --model ${params.model}" }
    if ( params.missing ) { extra_flags += " --missing ${params.missing}" }
    if ( params.no_regress ) { extra_flags += " --no-regress ${params.no_regress}" }
    if ( params.score ) { extra_flags += " --score ${params.score}" }
    if ( params.upper ) { extra_flags += " --upper ${params.upper}" }

    // 7 - R specific commands not available (pipeline is using a Docker image that handles such details)

    // 8 - Plotting

    if ( params.bar_col_high ) { extra_flags += " --bar-col-high ${params.bar_col_high}" }
    if ( params.bar_col_low ) { extra_flags += " --bar-col-low ${params.bar_col_low}" }
    if ( params.bar_col_p ) { extra_flags += " --bar-col-p ${params.bar_col_p}" }
    if ( params.bar_palatte ) { extra_flags += " --bar-palatte ${params.bar_palatte}" }
    if ( params.multi_plot ) { extra_flags += " --multi-plot ${params.multi_plot}" }
    if ( params.plot ) { extra_flags += " --plot ${params.plot}" }
    if ( params.plot_set ) { extra_flags += " --plot-set ${params.plot_set}" }
    if ( params.scatter_r2 ) { extra_flags += " --scatter-r2 ${params.scatter_r2}" }

    // 8 - Miscellaneous

    if ( params.all_score ) { extra_flags += " --all-score ${params.all_score}" }
    if ( params.exclude ) { extra_flags += " --exclude ${params.exclude}" }
    if ( params.extract ) { extra_flags += " --extract ${params.extract}" }
    if ( params.ignore_fid ) { extra_flags += " --ignore-fid ${params.ignore_fid}" }
    if ( params.keep_ambig ) { extra_flags += " --keep-ambig ${params.keep_ambig}" }
    if ( params.logit_perm ) { extra_flags += " --logit-perm ${params.logit_perm}" }
    if ( params.memory ) { extra_flags += " --memory ${params.memory}" }
    if ( params.non_cumulate ) { extra_flags += " --non-cumulate ${params.non_cumulate}" }
    if ( params.out ) { extra_flags += " --out ${params.out}" }
    if ( params.perm ) { extra_flags += " --perm ${params.perm}" }
    if ( params.print_snp ) { extra_flags += " --print-snp ${params.print_snp}" }
    if ( params.seed ) { extra_flags += " --seed ${params.seed}" }
    if ( params.thread ) { extra_flags += " --thread ${params.thread}" }
    if ( params.ultra ) { extra_flags += " --ultra ${params.ultra}" }
    if ( params.x_range ) { extra_flags += " --x-range ${params.x_range}" }
    
    quantile_flag = params.quantile =~ false ? '' : "--quantile ${params.quantile}"
    quant_break_flag = params.quant_break =~ false ? '' : "--quant-break ${params.quant_break}"
    quant_ref_flag = params.quant_ref =~ false ? '' : "--quant-ref ${params.quant_ref}"
    quant_extract_flag = params.quant_extract =~ false ? '' : "--quant-extract ${params.quant_extract}"
    bgen_flag = params.gt_format == "dosages" ? "--type bgen" :  ""
    '''
    if [ "!{bgen_flag}" != "" ]; then
        ls *.bgen > genotype_files.txt
    else
        ls *.bed > genotype_files.txt
    fi
    while IFS= read -r line; do echo \${line%%.*} >> target_list.txt; done < genotype_files.txt
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base !{base} \\
        --snp SNPID \\
        --chr CHR \\
        --bp POS \\
        --A1 Allele1 \\
        --A2 Allele2 \\
        --pvalue p.value \\
        --stat BETA \\
        !{bgen_flag} \\
        --beta \\
        --target-list target_list.txt \\
        --binary-target !{params.binary_trait} \\
        --pheno !{pheno} \\
        --cov !{cov} !{extra_flags} !{quantile_flag} !{quant_break_flag} !{quant_ref_flag} !{quant_extract_flag} \\

    # Remove date from image names (only for images produced by PRSice)
    images=$(ls *.png)
    for image in $images; do
        date=$(echo $image | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
        if ! [[ -z "${date// }" ]]; then
        mv "${image}" "${image/${date}/}"
        fi
    done

    # Remove date for quantile table (if quantile plot was produced)
    if ls PRSice_QUANTILES*.txt 1> /dev/null 2>&1; then
        table=$(ls PRSice_QUANTILES*.txt)
        date=$(echo $table | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
        if ! [[ -z "${date// }" ]]; then
        mv "${table}" "${table/${date}/}"
        fi
    fi
    '''
}

process additional_plots {

    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(pheno), path(cov)
    path(prs)
    path(metadata)

    output:
    path("*.png"), emit: more_plots

    script:
    """
    plot_prs_vs_cov.R --input_cov ${cov} --input_prs ${prs} --input_metadata ${metadata}
    plot_prs_density.R --input_pheno ${pheno} --input_prs ${prs}
    """
}

process produce_report {

    publishDir params.outdir, mode: "copy"

    input:
    path(plots)
    path("*")

    output:
    path("MultiQC/multiqc_report.html"), emit: reports

    script:
    if (params.quantile) {
        quantile_plot = "PRSice_QUANTILES_PLOT.png"
        quantile_table = "PRSice_QUANTILES.txt"
    } else {
        quantile_plot = "FALSE"
        quantile_table = "FALSE"
    }
    """
    cp /opt/bin/* .

    R -e "rmarkdown::render('prs_report.Rmd', params = list(barplot='PRSice_BARPLOT.png', highres.plot='PRSice_HIGH-RES_PLOT.png', density.plot='prs-density.png', quantile.plot='${quantile_plot}', quantile.table='${quantile_table}', prs.prsice='PRSice.prsice', prs.summary='PRSice.summary'))"
    mkdir MultiQC && mv prs_report.html MultiQC/multiqc_report.html

    """
}

process ldpred_coord {

    label "ldpred1"
    publishDir "${params.outdir}/LDpred", mode: "copy"
        
    input:
    each path(base)
    each path(merged_plink_file)
        
    output:
    path("ldpred.coord"), emit: harmonised_coord

    shell:
    '''
    sed "s/ /\\t/g" !{base}  > test_sumstats_tab.tsv
    ldpred coord \
    --chr CHR \
    --out ldpred.coord \
    --gf merged \
    --ssf-format CUSTOM \
    --ssf test_sumstats_tab.tsv \
    --pos POS \
    --A1 Allele1 \
    --A2 Allele2 \
    --pval p.value \
    --eff BETA \
    --se SE \
    --N 1000 \
    --rs SNPID 1> ldpred_coord.log
    '''
}

process ldpred_gibbs {

    label "ldpred1"
    publishDir "${params.outdir}/LDpred", mode: "copy"
        
    input:
    each path(cf_file)
    each path(merged_plink_file)
        
        
    output:
    path("ldpred.weights*"), emit: ldpred_weights

    script:
    """
    ldpred gibbs \
    --cf ${cf_file} \
    --ldr 150 \
    --f 1 0.3 0.1 0.03 0.01 \
    --out ldpred.weights \
    --ldf merged 1> ldpred.weights.log
    """
}

process ldpred_score {

    label "ldpred1"
    publishDir "${params.outdir}/LDpred", mode: "copy"
        
    input:
    each path(ldpred_weights)
    each path(merged_plink_file)
    tuple path(pheno), path(cov)
        
    output:
    path("ldpred.scores*"), emit: ldpred_scores
    path("ldpred.score_LDpred-inf.txt"), emit: scores_to_percentiles

    script:
    """
    awk '{print \$1,\$6}' merged.fam >> phenotypes.txt
    ldpred score \
    --gf merged \
    --rf ldpred.weights \
    --f 1 0.3 0.1 0.03 0.01 \
    --pf phenotypes.txt \
    --pf-format STANDARD \
    --out ldpred.score 1> ldpred.scores.log

    for scores in ldpred*.txt; do sed -i 's/ //g' \${scores}; done  # remove white spaces from scores

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


workflow lifebitai_prs_multi_plink_saige{
        take:
            ch_saige_base
            ch_target_pheno

        main:
            if (params.prsice) {
                /*------------------------------------------------------------------
                Setting up target pheno file
                ---------------------------------------------------------------------*/
                pheno_metadata_ch = Channel
                    .fromPath(params.pheno_metadata, checkIfExists: true)
                    .ifEmpty { exit 1, "Phenotype metadata file not found: ${params.pheno_metadata}" }

            }
            /*-----------------------------------------------------------------------------------------
            Setting up target dataset: PLINK files and pheno file output from lifebit-ai/biobank-gwas
            --------------------------------------------------------------------------------------------*/

            if (params.ldpred1 || params.prsice) {
                if (params.target_plink_dir) {
                    if (params.gt_format=="hard_called") {
                        ch_target_genotypic_dir = Channel
                            .fromPath("${params.target_plink_dir}/*.{bed,bim,fam}")
                            .ifEmpty { error "No target plink files found in : ${params.target_plink_dir}" }
                            .map { file ->
                                def key = file.name.toString().tokenize('_').get(0)
                                return tuple(key, file)
                            }
                            .groupTuple()

                    } else if (params.gt_format == "dosages") {
                        ch_target_genotypic_dir = Channel
                            .fromPath("${params.target_plink_dir}/*.{bgen,bgen.bgi,sample}")
                            .ifEmpty { error "No target bgen/bgen.bgi/.sample files found in : ${params.target_plink_dir}" }
                            .map { file ->
                                def key = file.name.toString().tokenize('_').get(0)
                                return tuple(key, file)
                            }
                            .groupTuple()
                    }

                }
            }

            ch_prs_scores_tables = Channel.empty()

            projectDir = workflow.projectDir


            /*----------------------------
            Setting up other parameters
            ------------------------------*/

            // // Initialise variable to store optional parameters
            transform_saige_base(ch_saige_base)

            /*---------------------------------
            Transforming target pheno input 
            -----------------------------------*/
            transform_target_pheno(ch_target_pheno)


            if (params.prsice) {
                polygen_risk_calcs( transform_saige_base.out.transformed_base,
                                    ch_target_genotypic_dir,
                                    transform_target_pheno.out.transformed_target_pheno )

                /*--------------------------
                Additional visualizations
                ----------------------------*/
                additional_plots( transform_target_pheno.out.transformed_target_pheno,
                                    polygen_risk_calcs.out.best_PRS,
                                    pheno_metadata_ch )

                /*--------------------------
                Produce R Markdown report                          
                ----------------------------*/
                produce_report( additional_plots.out.more_plots,
                                polygen_risk_calcs.out.results_for_report )

                ch_scores_to_percentiles = polygen_risk_calcs.out.best_PRS
            }

            if ( params.ldpred1 ) {
                merge_plink(ch_target_genotypic_dir)

                ldpred_coord(transform_saige_base.out.transformed_base_ldpred,
                                merge_plink.out.merged_plink)
                
                ldpred_gibbs(ldpred_coord.out.harmonised_coord,
                                merge_plink.out.merged_plink)
                
                ldpred_score(ldpred_gibbs.out.ldpred_weights,
                                merge_plink.out.merged_plink,
                                transform_target_pheno.out.transformed_target_pheno)
                
                ch_scores_to_percentiles = ldpred_score.out.scores_to_percentiles
            }

            if (params.calc_prs_percentiles) {
                calculate_prs_percentiles(ch_scores_to_percentiles)
            }

            if (params.filter_by_percentiles) {
                filter_by_percentiles_v1(merge_plink.out.merged_plink,
                                            calculate_prs_percentiles.out.samples_filtered_by_percentile)
            
                ch_filtered_percentile_cohort = filter_by_percentiles_v1.out.filtered_percentile_cohort
            }
        
        emit:
            ch_scores_to_percentiles
}

workflow{
    allowed_gt_list = allowed_genotypic_formats()
    // Check inputs

    ch_target_pheno = Channel
            .fromPath(params.target_pheno, checkIfExists: true)
            .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
            
    ch_saige_base = Channel
                .fromPath(params.saige_base, checkIfExists: true)
                .ifEmpty { exit 1, "SAIGE summary stats (base cohort) not found: ${params.saige_base}" }

     // Check that one and only one tool is called at a time
    tools_called = (params.ldpred1 ? 1 : 0) + \
                (params.prsice ? 1 : 0) 
    if (tools_called != 1) {
        exit 1, """Please select one, and only one, of the included tools to run:
                --ldpred1
                --prsice"""
    }

    lifebitai_prs_multi_plink_saige(ch_saige_base,
                                ch_target_pheno)

}

