#!/usr/bin/env nextflow
/*
========================================================================================
                         iudlgnt-gwas-sumstats-utils
========================================================================================
iudlgnt-gwas-sumstats-utils
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*-----------
  Processes  
--------------*/

// Do not delete this process
// Create introspection report

process obtain_pipeline_metadata {
  publishDir "${params.tracedir}", mode: "copy"

  input:
  val repository
  val commit
  val revision
  val script_name
  val script_file
  val project_dir
  val launch_dir
  val work_dir
  val user_name
  val command_line
  val config_files
  val profile
  val container
  val container_engine
  val raci_owner
  val domain_keywords

  output:
  path("pipeline_metadata_report.tsv"), emit: pipeline_metadata_report
  
  shell:
  '''
  echo "Repository\t!{repository}"                  > temp_report.tsv
  echo "Commit\t!{commit}"                         >> temp_report.tsv
  echo "Revision\t!{revision}"                     >> temp_report.tsv
  echo "Script name\t!{script_name}"               >> temp_report.tsv
  echo "Script file\t!{script_file}"               >> temp_report.tsv
  echo "Project directory\t!{project_dir}"         >> temp_report.tsv
  echo "Launch directory\t!{launch_dir}"           >> temp_report.tsv
  echo "Work directory\t!{work_dir}"               >> temp_report.tsv
  echo "User name\t!{user_name}"                   >> temp_report.tsv
  echo "Command line\t!{command_line}"             >> temp_report.tsv
  echo "Configuration file(s)\t!{config_files}"    >> temp_report.tsv
  echo "Profile\t!{profile}"                       >> temp_report.tsv
  echo "Container\t!{container}"                   >> temp_report.tsv
  echo "Container engine\t!{container_engine}"     >> temp_report.tsv
  echo "RACI owner\t!{raci_owner}"                 >> temp_report.tsv
  echo "Domain keywords\t!{domain_keywords}"       >> temp_report.tsv

  awk 'BEGIN{print "Metadata_variable\tValue"}{print}' OFS="\t" temp_report.tsv > pipeline_metadata_report.tsv
  '''
}

process fetch_from_ieu {
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/ieu/", mode: 'copy'
    }

    input:
    val(study_id)

    output:
    tuple val("${study_id}"), path("${study_id}.vcf.gz"), path("${study_id}_metadata_ieu.tsv"), emit: sumstats
    path("reported_traits.tsv"), emit: reported_traits

    script:
    """
    ieu_fetch.sh $study_id
    """
}

process munge_ieu {
    label 'munge'
    label "high_memory"
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/ieu/", mode: 'copy'
    }

    input:
    tuple val(study_id), path(sumstats), path(metadata)
    each path(snplocs_grch37)
    each path(snplocs_grch38)

    output:
    tuple val(study_id), path("${study_id}_harmonised_sumstats.vcf"), emit: harmonised_sumstats

    script:
    """
    R CMD INSTALL $snplocs_grch37
    R CMD INSTALL $snplocs_grch38

    munge_ieu.R \
        --gwas_vcf $sumstats \
        --meta_table $metadata \
        --method ${params.munge_method} \
        --study_id $study_id \
        --output "${study_id}_harmonised_sumstats.vcf" \
        --col_miss ${params.miss_percent_allow} \
        --dbsnp ${params.dbsnp} \
        --field_descriptions field_descriptions.tsv \
        --ncpus ${task.cpus}
    
    gunzip -S .bgz ${study_id}_harmonised_sumstats.vcf.bgz
    """
}


process standardisation_beta2or {
    label 'python'
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/conversion/", mode: 'copy'
    }

    input:
    tuple val(study_id), path(sumstats)

    output:
    tuple val("${study_id}"), path("${study_id}_harmonised_conv_sumstats.vcf"), emit: conv_sumstats
    
    script:
    """
    # Detect the FRQ column necessary to `--standardise`
    FRQ_COLUMN_EXISTS=\$(awk '/^[^#]/ {if(\$9 ~ /AF/) a=1; else a=0; exit}END{print a}' $sumstats)

    # Decide and convert
    if [ \$FRQ_COLUMN_EXISTS == "1" ]; then
        convert_coeff.py -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --standardise --beta2or
        echo "[MSG] BETA and SE standardisation performed."
        echo "[MSG] BETA to OR conversion performed."
    else
        convert_coeff.py -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --beta2or
        echo "[MSG] BETA to OR conversion performed."
        echo "[WARNING] BETA and SE standardisation not performed as FRQ value not present."
    fi
    """
}

process standardisation {
    label 'python'
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/conversion/", mode: 'copy'
    }

    input:
    tuple val(study_id), path(sumstats)

    output:
    tuple val("${study_id}"), path("${study_id}_harmonised_conv_sumstats.vcf"), emit: conv_sumstats
    
    script:
    """
    # Detect the FRQ column necessary to `--standardise`
    FRQ_COLUMN_EXISTS=\$(awk '/^[^#]/ {if(\$9 ~ /AF/) a=1; else a=0; exit}END{print a}' $sumstats)

    # Decide and convert
    if [ \$FRQ_COLUMN_EXISTS == "1" ]; then
        convert_coeff.py -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --standardise
        echo "[MSG] BETA and SE standardisation performed."
    else
        cp $sumstats ${study_id}_harmonised_conv_sumstats.vcf
        echo "[WARNING] BETA and SE standardisation not performed as FRQ value not present."
    fi
    """
}

process beta2or {
    label 'python'
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/conversion/", mode: 'copy'
    }

    input:
    tuple val(study_id), path(sumstats)

    output:
    tuple val("${study_id}"), path("${study_id}_harmonised_conv_sumstats.vcf"), emit: conv_sumstats
    
    script:
    """
    # BETA value should be always present
    convert_coeff.py -v $sumstats -o ${study_id}_harmonised_conv_sumstats.vcf --beta2or
    """
}

process map_traits {
    label "mid_memory"
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/metadata/", mode: 'copy'
    }

    input:
    path('reported_traits.list')
    path(omop_vocabulary)

    output:
    path("mapped.csv"), emit: mapped_traits

    script:
    """
    unzip ${omop_vocabulary}

    echo "traitName" > traits.csv
    cat reported_traits.list >> traits.csv

    map_traits.R --traits_csv=traits.csv --vocabulary_folder=omop-vocab-files --snomed_grouping_level=4
    """
}

process insert_traits {
    if (params.keep_intermediate_files) {
        publishDir "${params.outdir}/metadata/", mode: 'copy'
    }

    input:
    tuple val(study_id), path(sumstats_vcf)
    each path("mapped.csv")
    each path('field_descriptions.tsv')

    output:
    tuple val(study_id), path("${study_id}_harmonised_conv_omop_sumstats.vcf"), emit: omop_gwas_vcf

    script:
    """
    sed '/^#CHROM/ q' $sumstats_vcf > header.txt

    add_trait_metadata.R \
        --gwas_vcf_head header.txt \
        --trait_table mapped.csv \
        --field_descriptions field_descriptions.tsv

    sed '/^#/ d' $sumstats_vcf  \
    | cat updated_header.txt - \
    > ${study_id}_harmonised_conv_omop_sumstats.vcf
    """
}

process make_qc_plots {
    label "mid_memory"
    publishDir "${params.outdir}/QC_plots", mode: 'copy'
    
    input:
    tuple val(study_id), path("fully_harmonised_sumstats.vcf")

    output:
    tuple val(study_id), path("${study_id}/*"), emit: qc_plots

    script:
    """
    # Check whether the VCF contains FRQ column
    if `tail -1 fully_harmonised_sumstats.vcf | grep -q "AF"`; then
        gwas_qc_plots.py \
            --sumstat fully_harmonised_sumstats.vcf \
            --c-maf AF \
            --c-beta ES \
            --c-se SE \
            --c-pval LP \
            --out_dir ./ \
            --study_identifier ${study_id}
    else
        mkdir ${study_id}
        touch ${study_id}/${study_id}_no_plots.txt
    fi
    """
}

process filter_sumstats {
    label 'bcftools'
    publishDir "${params.outdir}/harmonised/filtered/", mode: 'copy'
    echo true
    
    input:
    tuple val(study_id), path(sumstats)

    output:
    tuple val(study_id), path("${study_id}_harmonised_filtered_sumstats.vcf"), emit: filtered_sumstats

    script:
    beta_smaller_filter = params.filter_beta_smaller_than || params.filter_beta_smaller_than == 0 ? "-e FORMAT/ES[*]<$params.filter_beta_smaller_than" : ""
    beta_greater_filter = params.filter_beta_greater_than || params.filter_beta_greater_than == 0 ? "-e FORMAT/ES[*]>$params.filter_beta_greater_than" : ""
    p_filter = params.filter_p_greater_than || params.filter_p_greater_than == 0 ? "-e FORMAT/LP[*]>$params.filter_p_greater_than" : ""
    frq_smaller_filter = params.filter_freq_smaller_than || params.filter_freq_smaller_than == 0 ? "-e FORMAT/AF[*]<$params.filter_freq_smaller_than" : ""
    frq_greater_filter = params.filter_freq_greater_than || params.filter_freq_greater_than == 0 ? "-e FORMAT/AF[*]>$params.filter_freq_greater_than" : ""
    info_filter = params.filter_info ? "-e FORMAT/SI[*]<$params.filter_info" : ""
    missing_info_filter = params.filter_missing_info || params.filter_missing_info == 0 ? "-e FORMAT/SI[*]='.'" : ""
    """
    cp $sumstats temp_sumstats.vcf
    # Filter expressions to pass into bcftools filter
    # Check for the AF and SI columns
    if tail -1 temp_sumstats.vcf | grep -q "AF"; then
        expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter" "$frq_smaller_filter" "$frq_greater_filter")
    else
        echo "[WARNING] AF column not found, all allele frequency filters will be ignored"
        expression_list=("$beta_smaller_filter" "$beta_greater_filter" "$p_filter")
    fi
    if tail -1 temp_sumstats.vcf | grep -q "SI"; then
        expression_list+=("$info_filter" "$missing_info_filter")
    else
        echo "[WARNING] SI column not found, all imputation accuracy filters will be ignored"
    fi

    # Apply filters sequentially, otherwise bcftools would apply OR logic
    for e in "\${expression_list[@]}"; do
        if [ "\$e" != "" ]; then
            bcftools filter "\$e" temp_sumstats.vcf > temp_sumstats_new.vcf
            mv temp_sumstats_new.vcf temp_sumstats.vcf
        fi
    done
    mv temp_sumstats.vcf "${study_id}_harmonised_filtered_sumstats.vcf"
    """
}

process collapse_vcf {
    label 'python'
    publishDir "${params.outdir}/harmonised/", mode: 'copy'

    input:
    tuple val(study_id), path(gwas_vcf)

    output:
    path("${study_id}.harmonised.gwas.vcf")
    
    script:
    """
    collapse_vcf.py $gwas_vcf -o ${study_id}.harmonised.gwas.vcf
    """
}

process vcf2hail {
    publishDir "${params.outdir}/harmonised/hail/", mode: 'copy'

    input:
    tuple val(study_id), path(sumstats)

    output:
    tuple val(study_id), path("${study_id}_hail_matrix.mt/*"), emit: hail_sumstats

    script:
    """
    GENOME=\$( grep -e '##GenomeBuild=' $sumstats | cut -d "=" -f 2 )
    vcf2hail.py -v $sumstats -g \$GENOME -o "${study_id}_hail_matrix.mt"
    """
}

process create_report {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'
  
  input:
  tuple val(study_id), path("${study_id}/*")
  path(report_dir)

  output:
  path("${study_id}_multiqc_report.html"), emit: report_outputs_all
  path("multiqc_report.html"), emit: report_outputs1

  script:
  """
  cp -r ${report_dir}/* .
  cp ${study_id}/* .
  # Skip when no plots were generated due to lack of FRQ column
  if [ -f "${study_id}_no_plots.txt" ]; then
    touch multiqc_report.html
    touch ${study_id}_multiqc_report.html
  else
    for i in `ls *.png`; do name=`basename \$i .png | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
    for i in `ls *.tsv`; do name=`basename \$i .tsv | sed 's/-/_/g'`; echo "\$name='\$i'" >> file_list.txt;done
    cat file_list.txt | tr "\n" "," | sed 's/,\$//g' > file_list1.txt
    # Generates the report
    Rscript -e "rmarkdown::render('report.Rmd', params = list(`cat file_list1.txt`))"
    cp report.html multiqc_report.html
    mv report.html ${study_id}_multiqc_report.html
  fi
  """
}

workflow lifebitai_harmonisation_ieu{
    take:
        ch_input

    main:
        /*--------------------------------------------------------
        Defining and showing header with all params information 
        ----------------------------------------------------------*/

        // Header log info

        def summary = [:]

        if (workflow.revision) summary['Pipeline Release'] = workflow.revision

        summary['Output dir']                                  = params.outdir
        summary['Launch dir']                                  = workflow.launchDir
        summary['Working dir']                                 = workflow.workDir
        summary['Script dir']                                  = workflow.projectDir
        summary['User']                                        = workflow.userName
        summary['Standardise BETA and SE']                     = params.standardise
        summary['Coefficient conversion']                      = params.coef_conversion
        summary['EBI missingness percent allowed']             = params.miss_percent_allow
        summary['Keep intermediate files']                     = params.keep_intermediate_files
        summary['OMOP vocabulary DB']                          = params.omop_vocabulary
        summary['Convert to Hail']                             = params.convert_to_hail
        summary['Filter BETA smaller than']                    = params.filter_beta_smaller_than
        summary['Filter BETA greater than']                    = params.filter_beta_greater_than
        summary['Filter P value greater than']                 = params.filter_p_greater_than
        summary['Filter Alt AF smaller than']                  = params.filter_freq_smaller_than
        summary['Filter Alt AF greater than']                  = params.filter_freq_greater_than
        summary['Filter missing SI']                           = params.filter_missing_info
        summary['Filter INFO/SI smaller than']                 = params.filter_info

        log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
        log.info "-\033[2m--------------------------------------------------\033[0m-"



        /*-------------------------------------------------
        Setting up introspection variables and channels  
        ----------------------------------------------------*/

        // Importantly, in order to successfully introspect:
        // - This needs to be done first `main.nf`, before any (non-head) nodes are launched. 
        // - All variables to be put into channels in order for them to be available later in `main.nf`.

        ch_repository         = Channel.of(workflow.manifest.homePage)
        ch_commitId           = Channel.of(workflow.commitId ?: "Not available is this execution mode. Please run 'nextflow run ${workflow.manifest.homePage} [...]' instead of 'nextflow run main.nf [...]'")
        ch_revision           = Channel.of(workflow.manifest.version)

        ch_scriptName         = Channel.of(workflow.scriptName)
        ch_scriptFile         = Channel.of(workflow.scriptFile)
        ch_projectDir         = Channel.of(workflow.projectDir)
        ch_launchDir          = Channel.of(workflow.launchDir)
        ch_workDir            = Channel.of(workflow.workDir)
        ch_userName           = Channel.of(workflow.userName)
        ch_commandLine        = Channel.of(workflow.commandLine)
        ch_configFiles        = Channel.of(workflow.configFiles)
        ch_profile            = Channel.of(workflow.profile)
        ch_container          = Channel.of(workflow.container)
        ch_containerEngine    = Channel.of(workflow.containerEngine)

        /*----------------------------------------------------------------
        Setting up additional variables used for documentation purposes  
        -------------------------------------------------------------------*/
        ch_raci_owner      = Channel.of(params.raci_owner)
        ch_domain_keywords = Channel.of(params.domain_keywords)

        /*----------------------
        Setting up input data  
        -------------------------*/
        // Check if any filter was selected
        if (params.filter_beta_smaller_than || params.filter_beta_greater_than || params.filter_p_greater_than || params.filter_freq_smaller_than || params.filter_freq_greater_than || params.filter_missing_info || params.filter_info) {
            filter_activated = true
        } else {
            filter_activated = false
        }

        // Define channels from repository files

        projectDir = workflow.projectDir

        // Channel for omop vocabulary
        ch_omop_vocabulary = Channel.value(file(params.omop_vocabulary))

        // Channels for SNPlocs
        ch_snplocs_grch38 =  Channel.value(file(params.snplocs_grch38))
        ch_snplocs_grch37 =  Channel.value(file(params.snplocs_grch37))

        // Channels for scripts
        fetch_from_ieu(ch_input, 
                        ch_ieu_script)
        
        munge_ieu(fetch_from_ieu.out.sumstats,
                    ch_munge_ieu_script,
                    ch_munge_funcs_script,
                    ch_harmonise_metadata_funcs_script,
                    ch_field_descriptions,
                    ch_snplocs_grch37,
                    ch_snplocs_grch38)
        
        ch_harmonised_sumstats = munge_ieu.out.harmonised_sumstats

        if (params.standardise && params.coef_conversion) {
            standardisation_beta2or(ch_harmonised_sumstats)
            
            ch_conv_sumstats = standardisation_beta2or.out.conv_sumstats

        } else if (params.standardise) {
            standardisation(ch_harmonised_sumstats)

            ch_conv_sumstats = standardisation.out.conv_sumstats

        } else if (params.coef_conversion) {
            beta2or(ch_harmonised_sumstats)

            ch_conv_sumstats = beta2or.out.conv_sumstats

        } else {
            ch_conv_sumstats = ch_harmonised_sumstats
        }  

        if (params.map_traits) {
            // combine all traits into a single csv for improved efficiency
            //(loading omop vocabulary is the majority of the processing here)

            ch_collected_traits = ch_reported_traits
                .map{it.text}
                .collect()
                .map{it.join('\n')}

            map_traits(ch_collected_traits,
                        ch_omop_vocabulary)

            
            insert_traits(ch_conv_sumstats,
                            map_traits.out.mapped_traits,
                            ch_field_descriptions)
            
            ch_omop_gwas_vcf = insert_traits.out.omop_gwas_vcf
    
        } else {
            ch_omop_gwas_vcf = ch_conv_sumstats
        }

        make_qc_plots(ch_omop_gwas_vcf)

        if (filter_activated) {
            filter_sumstats(ch_omop_gwas_vcf)
            
            ch_filtered_sumstats = filter_sumstats.out.filtered_sumstats

        } else {
            ch_filtered_sumstats = ch_omop_gwas_vcf
        }

        collapse_vcf(ch_filtered_sumstats)

        if (params.convert_to_hail) {
            vcf2hail(ch_filtered_sumstats)
        }

        ch_report_dir = Channel.value(file("${projectDir}/bin/report"))

        create_report(make_qc_plots.out.qc_plots,
                        ch_report_dir)
        
        report_output = create_report.out.report_outputs_all

    emit:
        ch_conv_sumstats
        report_output
    }

workflow {
    // Define Channels from input
    ch_input = Channel.fromPath(params.input)
        .ifEmpty { exit 1, "Cannot find input file containing study IDs: ${params.input}" }
        .splitCsv()
        .flatten()
        .take(params.take_n_studies) //default is -1 i.e. take all files (but param is useful for testing with fewer files)

    harmonisation_ieu(ch_input)
}
