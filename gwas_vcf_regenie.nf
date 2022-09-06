#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/gwas
========================================================================================
 lifebit-ai/gwas GWAS pipeline using SAIGE linear mixed model approach for association testing
 #### Homepage / Documentation
 https://github.com/lifebit-ai/gwas
----------------------------------------------------------------------------------------
*/


/*---------------------------------------------------
  Define and show header with all params information
-----------------------------------------------------*/
nextflow.enable.dsl=2
// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

//summary['Max Resources']                  = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']                     = params.outdir
summary['Launch dir']                     = workflow.launchDir
summary['Working dir']                    = workflow.workDir
summary['Script dir']                     = workflow.projectDir
summary['User']                           = workflow.userName
summary['phenotype_colname']              = params.phenotype_colname
summary['pheno_transform']                = params.pheno_transform

summary['input_folder_location']          = params.input_folder_location
summary['file_pattern']                   = params.file_pattern
summary['file_suffix']                    = params.file_suffix
summary['index_suffix']                   = params.index_suffix
summary['number_of_files_to_process']     = params.number_of_files_to_process

summary['q_filter']                       = params.q_filter
summary['miss_test_p_threshold']          = params.miss_test_p_threshold
summary['variant_missingness']            = params.miss
summary['variant_minor_allele_freq']      = params.maf
summary['variant_minor_allele_count']     = params.mac
summary['hwe_threshold']                  = params.hwe_threshold
summary['extracted pruned region']        = params.extracted_prune_region

summary['remove_related_samples']         = params.remove_related_samples
summary['king_coefficient']               = params.king_coefficient
summary['king_plink_memory']              = params.king_plink_memory

summary['king_reference_data']            = params.king_reference_data
summary['run_ancestry_inference']         = params.run_ancestry_inference
summary['min_subpop_size']                = params.min_subpop_size

summary['run_pca']                        = params.run_pca
summary['number_pcs']                     = params.number_pcs
summary['remove_outliers_maxiter']           = params.remove_outliers_maxiter
summary['remove_outliers_sigma']             = params.remove_outliers_sigma


summary['gwas_cat']                       = params.gwas_cat

summary['outdir']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



/*--------------------------------------------------
Pre-GWAS filtering - download, filter and convert VCFs
---------------------------------------------------*/
process standardise_phenofile_and_get_samples {
  input:
  file('original.pheno.tsv')

  output:
  path("notransform.phe"), emit: pheno_no_transform
  path("all_samples.tsv"), emit: samples_vcf2plink

  script:
  """
  # make dummy transform file to move params.phenotype_colname to 3rd column
  # and IGNORE anything that is not the phenotype column or specified covariate_column

  if [ "${params.phenotype_colname}" = "false" ]; then
    pheno_col=\$(head -n 1 original.pheno.tsv | cut -f3 )
  else
    pheno_col=${params.phenotype_colname}
  fi

  if [ "${params.covariate_cols}" = "ALL" ]; then
    covar_cols=\$(head -n 1 original.pheno.tsv | cut --complement -f1,2 | tr '\\t' ',')
  elif [ "${params.covariate_cols}" = "NONE" ]; then
    covar_cols=" "
  else
    covar_cols=${params.covariate_cols}
  fi

  awk -v covar_cols="\$covar_cols" -v pheno_col="\$pheno_col" '\
    BEGIN{FS=OFS="\\t"} \
    NR==1{ \
      \$1="run_id"; \$2="test"; \
      print \$0; \
      split(covar_cols,covs,","); for(i in covs){cols[covs[i]]}; \
      cols[pheno_col]; \
      \$1="notransform"; \$2=pheno_col; \
      for (i=3; i <= NF; i++){if(\$i in cols){\$i=""}else{\$i="IGNORE"}}; \
      print \$0 \
    }' original.pheno.tsv > dummmy_transform.tsv

  transform_pheno.R \
    --pheno original.pheno.tsv \
    --transform dummmy_transform.tsv \
    --out_prefix ./

  cut -f1,2 notransform.phe > all_samples.tsv
  """
}

process vcf2plink {
  tag "$name"
  label 'high_memory'
  publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

  input:
  tuple val(name), val(chr), file(vcf), file(index)
  each file(phe_file)

  output:
  tuple val(name), val(chr), path('*.bed'), path('*.bim'), path('*.fam'), emit: filteredPlink

  script:
  plink_memory = extractInt(task.memory.toString()) * 1000
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  tail -n +2 ${phe_file}| cut -f2 > samples.txt
  bcftools view -S samples.txt $vcf -Oz -o ${name}_downsampled.vcf.gz
  bcftools view -q ${params.q_filter} ${name}_downsampled.vcf.gz -Oz -o ${name}_filtered.vcf.gz
  bcftools index ${name}_filtered.vcf.gz

  # Create PLINK binary from vcf.gz
  plink2 \
    --make-bed \
    --set-missing-var-ids @:#,\\\$r,\\\$a \
    --vcf ${name}_filtered.vcf.gz \
    --out ${name}_filtered \
    --vcf-half-call m \
    --memory ${plink_memory} \
    --double-id \
    --set-hh-missing \
    --new-id-max-allele-len 60 missing
  """
}

process filter_genotypic_data {
  tag "$name"
  publishDir "${params.outdir}/filter_genotypic_data", mode: 'copy'
  input:
  tuple val(name), val(chr), file(bed), file(bim), file(fam)

  output:
  tuple path("${name}.misHWEfiltered.bed"), path("${name}.misHWEfiltered.bim"), path("${name}.misHWEfiltered.fam"), emit: filtered_geno_data

  script:
  plink_memory = extractInt(task.memory.toString()) * 1000
  """
  ### PROCESS: basic_qa_filtering
  #### Basic pre-gwas QC:
  #### Filter out all variants with missingness > 0.1 [--geno 0.1]
  #### Filter out all variants with minor allele frequency < 0.05 [--maf 0.05]
  #### Filter out all variants with minor allele count < 1 [--mac 1]
  #### For reference: minor allele count - number of times minor allele was observed in the cohort. E.g. mac=237 would mean 237 participants had this allele.

  plink \
    --bfile ${name}_filtered \
    --allow-no-sex \
    --mind ${params.mind_threshold} \
    --geno ${params.miss} \
    --maf ${params.maf} \
    --mac ${params.mac} \
    --memory ${plink_memory} \
    --1 \
    --keep-allele-order \
    --make-bed \
    --out ${name}_QA_filtered

  ### PROCESS: filter_hwe
  # plink -> hwe_filtered plink
  plink \
    --bfile ${name}_QA_filtered \
    --memory ${plink_memory} \
    --allow-no-sex \
    --hwe ${params.hwe_threshold} ${params.hwe_test} \
    --out ${name}.misHWEfiltered \
    --make-bed \
    --1 \
    --keep-allele-order
  """
}

process merge_plink {
  publishDir "${params.outdir}/merged_plink", mode: 'copy'

  input:
  path("*")

  output:
  tuple path('merged.bed'), path('merged.bim'), path('merged.fam'), emit: plink_merged

  script:
  plink_memory = extractInt(task.memory.toString()) * 1000
  """
  ls *.bed > bed.txt
  ls *.bim > bim.txt
  ls *.fam > fam.txt
  paste bed.txt bim.txt fam.txt > merge.temp.list
  tail -n +2 merge.temp.list > merge.list
  bed_file=\$(head -n1 merge.temp.list | cut -f1)
  bed_prefix=`echo "\${bed_file%.*}"`
  plink --keep-allele-order \
  --bfile \${bed_prefix} \
  --merge-list merge.list \
  --allow-no-sex \
  --memory ${plink_memory} \
  --make-bed \
  --out merged

  """
}

process check_sample_sex {
  tag "Sex check"
  publishDir "${params.outdir}/sex_check", mode: 'copy'

  input:
  tuple file('merged.bed'), file('merged.bim'), file('merged.fam')

  output:
  path("sex_check.log"), emit: sex_check_log

  script:
  """
  x_chrom=\$(cut -f1 merged.bim | grep 23 || true )

  if [ -z \$x_chrom ]; then
    echo "No chromosome X detected in genotype data - continuing without performing sex check/sex imputation." > sex_check.log
  else
    echo "X chromosome identified in genotype data - performing sex check." > sex_check.log
    plink --keep-allele-order --bfile merged --check-sex --out sex_check
  fi

  """
  }

process ld_prune {
  publishDir "${params.outdir}/merged_pruned_plink", mode: 'copy'

  input:
  tuple file('merged.bed'), file('merged.bim'), file('merged.fam')
  file(long_range_ld_regions)

  output:
  tuple val('merged_pruned'), path('merged_pruned.bed'), path('merged_pruned.bim'), path('merged_pruned.fam'), emit: pruned_variants //(ch_pruned_variants_for_grm, ch_pruned_variants_for_relatedness, ch_unrelated_for_pca, ch_unrelated_for_ancestry_inference)

  script:
  plink_memory = extractInt(task.memory.toString()) * 1000
  extract_options = params.extracted_prune_region ? "--extract merged.prune.in" : ""
  """
  plink \
  --bfile merged \
  --keep-allele-order \
  --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
  --exclude range ${long_range_ld_regions} \
  --allow-no-sex \
  --memory ${plink_memory} \
  --out merged

  plink \
  --bfile merged \
  --keep-allele-order \
  ${extract_options} \
  --make-bed \
  --memory ${plink_memory} \
  --allow-no-sex \
  --out merged_pruned
  """
}

process calculate_relatedness {
    label "tiny_memory"

    input:
    tuple val(name), file(bim), file(bed), file(fam)

    output:
    path("relatedness.king.cutoff.in.id"), emit: related_filter_keep_files

    script:
    plink_memory = params.king_plink_memory ? "--memory ${params.king_plink_memory}" : ""
    """
    plink2 \
        --bfile ${bim.baseName} \
        --king-cutoff ${params.king_coefficient} \
        $plink_memory \
        --out relatedness
    """
}

process infer_ancestry {
    label "tiny_memory"
    publishDir "${params.outdir}/king_ancestry_inference/", mode: 'copy'

    input:
    tuple file('keep.tsv'), val(name), file('in.bed'), file('in.bim'), file('in.fam')
    file(ancestry_inference_Rscript)
    set val(ref_name), file('ref.bed.xz'), file('ref.bim.xz'), file('ref.fam.xz')

    output:
    path("*.keep.tsv"), emit: split_ancestry_keep_files_out
    tuple path("*_InferredAncestry.txt"), path("*ancestryplot.*"), path("*pc.txt")

    script:
    """
    # decompress .xz suffixed king reference files
    unxz --force ref.bed.xz
    unxz --force ref.bim.xz
    unxz --force ref.fam.xz

    # subset input to just samples to keep
    cut -f 1-2 keep.tsv > inputkeep.tsv
    plink2 --bfile in --keep inputkeep.tsv --make-bed --out subset

    # Swap placeholders with user provided values
    sed -i "s/PREFIX_PLACEHOLDER/subset/g" Ancestry_Inference.R

    # calculate PCs
    king -b ref.bed,subset.bed --pca --projection --rplot --prefix subset --cpus ${task.cpus}

    # Run ancestry inference script
    Rscript Ancestry_Inference.R subsetpc.txt subset_popref.txt subset

    awk 'BEGIN{OFS="\\t"} NR>1{print \$1,\$2 > \$5 ".keep.tsv" }' subset_InferredAncestry.txt
    """
}

process filter_pca {
  tag "${ancestry_group}"
  publishDir "${params.outdir}/${ancestry_group}/pca/", mode: 'copy'

  input:
  tuple val(ancestry_group), file('in.tsv'), val(name), file('in.bed'), file('in.bim'), file('in.fam')

  output:
  tuple val(ancestry_group), path('pca_results_final.eigenvec'), path('pca_results_final.eigenval'), emit: pca_results
  tuple val(ancestry_group), path('out.keep.tsv'), emit: pca_keep_files
  path("remove_outliers_*.log"), emit: pca_logs

  script:
  """
  # subset to only the ancestry group we are currently operating on
  cut -f 1-2 in.tsv > in.keep.tsv
  plink2 --bfile in --keep in.keep.tsv --make-bed --out subset

  if [ \$(wc -l subset.bim | cut -d " " -f1) -lt 220 ]; then
      echo "Error: PCA requires data to contain at least 220 variants." 1>&2
      exit 1
  fi

  touch plink_remove_ids_init.tsv
  i=0
  n_outliers=-1

  while [ \$i -le ${params.remove_outliers_maxiter} ] && [ \$n_outliers -ne 0 ]; do

    cat plink_remove_ids_*.tsv > all_plink_remove_ids.tsv

    if [ \$(wc -l in.fam | cut -d " " -f1) -gt 5000 ]; then
      plink2 --bfile subset --remove all_plink_remove_ids.tsv --pca ${params.number_pcs} approx --out pca_results_\${i} 1>&2
    else
      plink2 --bfile subset --remove all_plink_remove_ids.tsv --pca ${params.number_pcs} --out pca_results_\${i} 1>&2
    fi

    pca_outliers.R --input=pca_results_\${i}.eigenvec \
      --out=plink_remove_ids_\${i}.tsv \
      --sigma=${params.remove_outliers_sigma} > remove_outliers_\${i}.log

    n_outliers=\$(wc -l plink_remove_ids_\${i}.tsv | cut -d " " -f1)
    i=\$[\$i+1]

  done

  i=\$[\$i-1]
  cp pca_results_\${i}.eigenval pca_results_final.eigenval
  cp pca_results_\${i}.eigenvec pca_results_final.eigenvec

  tail -n +2 pca_results_final.eigenvec | cut -f 1-2 > out.keep.tsv
  """
}

process transform_phenofile {
  input:
  file(original_pheno)
  file(transform_tsv)

  output:
  file("*.phe"), emit: transform_pheno_out

  script:
  """
  transform_pheno.R --pheno ${original_pheno} --transform ${transform_tsv} --out_prefix ./
  """

}

process create_ancestry_x_transform_pheno {
  tag "${ancestry_group} ${gwas_tag}"
  
  input:
  tuple val(gwas_tag), file('in.phenofile.phe')
  tuple val(ancestry_group), file('keep.tsv')

  output:
  tuple val(ancestry_group), val(gwas_tag), file('out.phenofile.phe'), emit: crossed_pheno_out

  script:
  """
  awk '\
    BEGIN{FS=OFS="\\t"} \
    ARGIND==1{a[\$1 "%" \$2]++; next} \
    FNR==1{print \$0; next} \
    a[\$1 "%" \$2]>0{print \$0} \
  ' keep.tsv in.phenofile.phe > out.phenofile.phe
  """
}

process add_pcs_to_pheno {
  tag "${ancestry_group} ${gwas_tag}"

  input:
  tuple val(ancestry_group), file('pca_results.eigenvec'), file('pca_results.eigenval'), val(gwas_tag), file('in.phenofile.phe')

  output:
  tuple val(ancestry_group), val(gwas_tag), path('out.phenofile.phe'), emit: final_phenos

  script:
  """
  awk '
    BEGIN{FS=OFS="\\t"} \
    ARGIND==1{fid=\$1;iid=\$2;\$1="%%";\$2="%%";a[fid "%" iid] = \$0} \
    ARGIND==2{print \$0, a[\$1 "%" \$2]} \
  ' pca_results.eigenvec in.phenofile.phe \
  | sed -r 's/%%\\t//g' > out.phenofile.phe
  """
}

process align_pheno_with_test_variant_plink {
  tag "${ancestry_group} ${gwas_tag}"
  publishDir "${params.outdir}/${ancestry_group}/${gwas_tag}/", mode: 'copy'

  input:
  tuple val(ancestry_group), val(gwas_tag), file('phenofile.phe')
  tuple file('in.bed'), file('in.bim'), file('in.fam')

  output:
  tuple val(ancestry_group), val(gwas_tag), file('trait_type'), file('phenofile.phe'), file('aligned_test_variants.bed'), file('aligned_test_variants.bim'), file('aligned_test_variants.fam'), emit: align_pheno_with_test_variant_plink_out

  script:
  """
  cut -f1-2 phenofile.phe > keep_samples.txt
  plink2 --bfile in --keep keep_samples.txt --make-bed --out aligned_test_variants
  
  uniquevals=\$(cut -f 3 phenofile.phe | tail -n +2 | sort -u | wc -l)
  if [ \$uniquevals -gt 2 ]; then
    echo -n "quantitative" > trait_type
  else
    echo -n "binary" > trait_type
  fi
  """
}

process filter_binary_missingness {
    tag "${ancestry_group} ${gwas_tag}"
    publishDir "${params.outdir}/${ancestry_group}/${gwas_tag}/", mode: 'copy'

    input:
    tuple val(ancestry_group), val(gwas_tag), val(trait_type), file('phenofile.phe'), file('in.bed'), file('in.bim'), file('in.fam')

    output:
    tuple val(ancestry_group), val(gwas_tag), val(trait_type), file('phenofile.phe'), file('*_filtered.bed'), file('*_filtered.bim'), file('*_filtered.fam'), emit: aligned_test_vars_binary_hwe_filtered

    script:
    plink_memory = extractInt(task.memory.toString()) * 1000
    """
    plink \
      --bfile in \
      --pheno phenofile.phe \
      --allow-no-sex \
      --test-missing midp \
      --memory ${plink_memory} \
      --out hwe \
      --keep-allele-order \

    awk '\$5 < ${params.miss_test_p_threshold} {print \$2 }' hwe.missing > hwe.missing_FAIL

    plink --bfile in \
      --keep-allele-order \
      --allow-no-sex \
      --exclude hwe.missing_FAIL \
      --memory ${plink_memory} \
      --make-bed \
      --out aligned_test_variants.binary_hwe_miss_filtered
    """
}

process convert2bgen   {
  tag "${ancestry_group} ${gwas_tag}"
  label "convert2bgen"

  input:
  tuple val(ancestry_group), val(gwas_tag), val(trait_type), file('phenofile.phe'), file('in.bed'), file('in.bim'), file('in.fam')

  output:
  tuple val(ancestry_group), val(gwas_tag), val(trait_type), path('phenofile.phe'), path('aligned_filtered.bgen'), path('aligned_filtered.sample'), emit: merged_bgen

  script:
  """
  plink2 --bfile in --export bgen-1.2 bits=8 --out aligned_filtered
  """
}

process regenie_step1_fit_model {
  tag "${ancestry_group} ${gwas_tag}"
  label 'regenie'
  publishDir "${params.outdir}/${ancestry_group}/${gwas_tag}/regenie/", mode: 'copy', pattern: "${ancestry_group}-${gwas_tag}-regenie_step1*"

  input:
  tuple val(ancestry_group), val(gwas_tag), val(trait_type), file('phenofile.phe'), file('in.bgen'), file('in.sample')

  output:
  tuple val(ancestry_group), val(gwas_tag), val(trait_type), path('phenofile.phe'), path('in.bgen'), path('in.sample'), path("*.loco"), path("*_pred.list"), path("covariates.txt"), path("pheno.txt"), emit: inputs_for_regenie_step2
  path("${ancestry_group}-${gwas_tag}-regenie_step1*")

  script:
  """
  sed -e '1s/^.//' phenofile.phe | sed 's/\t/ /g' > full_pheno_covariates.txt
  cut -d' ' -f1-3 full_pheno_covariates.txt > pheno.txt
  cut -d' ' --complement -f 3 full_pheno_covariates.txt > covariates.txt

  regenie \
  --step 1 \
    --bgen in.bgen \
    --covarFile covariates.txt \
    --phenoFile pheno.txt \
    --cc12 \
    --bsize 100 \
    --threads ${task.cpus} \
    ${trait_type == "binary" ? '--bt' : ''} \
    --lowmem \
    --lowmem-prefix tmp_rg \
    --use-relative-path \
    --out ${ancestry_group}-${gwas_tag}-regenie_step1
  """
}

process regenie_step2_association_testing {
  tag "${ancestry_group} ${gwas_tag}"
  label 'regenie'
  publishDir "${params.outdir}/${ancestry_group}/${gwas_tag}/regenie", mode: 'copy'

  input:
  tuple val(ancestry_group), val(gwas_tag), val(trait_type), path('phenofile.phe'), path('in.bgen'), path('in.sample'), file(loco), file(pred), file ("covariates.txt"), file ("pheno.txt")

  output:
  tuple val(ancestry_group), val(gwas_tag), file("${ancestry_group}-${gwas_tag}-regenie_firth*.regenie"), emit: regenie_out
  path("${ancestry_group}-${gwas_tag}-regenie_firth*"), emit: regenie_step2_assoc

  script:
  """
  regenie \
    --step 2 \
    --bgen in.bgen \
    --covarFile covariates.txt \
    --phenoFile pheno.txt \
    --cc12 \
    --bsize 200 \
    --threads ${task.cpus} \
    ${trait_type == "binary" ? '--bt' : ''} \
    --minMAC ${params.regenie_min_mac} \
    --minINFO ${params.regenie_min_imputation_score} \
    --firth --approx \
    --pThresh 0.01 \
    --pred ${pred} \
    --out ${ancestry_group}-${gwas_tag}-regenie_firth
  """
}

process munge_regenie {
  tag "${ancestry_group} ${gwas_tag}"
  label 'mungesumstats'
  publishDir "${params.outdir}/${ancestry_group}/${gwas_tag}/regenie", mode: 'copy'

  input:
  tuple val(ancestry_group), val(gwas_tag), file(regenie_table)

  output:
  tuple val(ancestry_group), val(gwas_tag), val('regenie'), path("${regenie_table.baseName}.vcf"), emit: summ_stats
  file("${regenie_table.baseName}.munge.log")

  script:
  """
  # rename GENPOS -> BP, remove last column, add col "P" calculated from col "LOG10P" 
  sed '/^##/d' $regenie_table \
  | sed '1 s/GENPOS/BP/' \
  | awk '{NF-=1} NR==1{print \$0, "P"; next} {print \$0, 10^(-\$13)}' \
  > regenie.tmp
  mungesumstats.R \
    --gwas_table=regenie.tmp \
    --outfile=${regenie_table.baseName}.vcf \
    ${params.map_pos2rsid ? "--map_pos=TRUE" : ""} \
    ${params.genome_build ? "--build=${params.genome_build}" : ""} \
    --ncpus=${task.cpus} 2> ${regenie_table.baseName}.munge.log
  """
}

/*--------------------------------------------------
  Defining functions
---------------------------------------------------*/

def extractInt( String input ) {
  return input.replaceAll("[^0-9]", "").toInteger()
}

def get_chromosome( file ) {
    // using RegEx to extract chromosome number from file name
    regexpPE = /(?:chr)[a-zA-Z0-9]+/
    (file =~ regexpPE)[0].replaceAll('chr','')
}

def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

def defineFormatList() {
    return [
      'bgen',
      'vcf'
    ]
}

workflow lifebitai_gwas_vcf_regenie{
  take:
    ch_user_input_vcf
    ch_king_reference_data
    ch_input_pheno_transform
    ch_high_ld_regions
    ch_gwas_cat
    ch_ld_scores
    ch_pheno

  main:
    standardise_phenofile_and_get_samples(ch_pheno)

    vcf2plink(ch_user_input_vcf,
              standardise_phenofile_and_get_samples.out.samples_vcf2plink)

    filter_genotypic_data(vcf2plink.out.filteredPlink)

    merge_plink(filter_genotypic_data.out.filtered_geno_data.collect())
    merged_plink = tuple "merged", merge_plink.out.plink_merged
    
    if (params.sex_check) {check_sample_sex(merge_plink.out.plink_merged)}

    ld_prune(merge_plink.out.plink_merged,
            ch_high_ld_regions)

    if (params.remove_related_samples) {
      calculate_relatedness(ld_prune.out.pruned_variants)
      ch_related_filter_keep_files = calculate_relatedness.out.related_filter_keep_files
    } 
    
    if (!params.remove_related_samples) {
      ch_related_filter_keep_files = ch_samples_for_no_related_filter
    }

    if (params.run_ancestry_inference) {
      ch_input_for_ancestry_inference = ch_related_filter_keep_files.combine(ch_unrelated_for_ancestry_inference)

      infer_ancestry(ch_input_for_ancestry_inference,
                      ch_king_reference_data)

      ch_ancestry_keep_files = infer_ancestry.out.split_ancestry_keep_files_out
                                  .flatMap()
                                  .map{[it.simpleName, it] }
                                  .filter{ it[1].readLines().size() >= params.min_subpop_size }
    } 

    if (!params.run_ancestry_inference) {
        ch_ancestry_keep_files = ch_related_filter_keep_files
          .map{["allancs"] + [it] }
    }

    if (params.run_pca) {
      ch_input_for_pca = ch_ancestry_keep_files.combine(ld_prune.out.pruned_variants)
      filter_pca(ch_input_for_pca)
      ch_pca_keep_files = filter_pca.out.pca_keep_files
    }

    if (!params.run_pca) {
      ch_pca_keep_files = ch_ancestry_keep_files
    }

    if (params.pheno_transform) {
      transform_phenofile(ch_pheno,
                          ch_input_pheno_transform)
      ch_transformed_phenos = transform_phenofile.out.transform_pheno_out
        .flatMap()
        .map{[it.baseName, it]}
    } 

    if (!params.pheno_transform) {
      ch_transformed_phenos = standardise_phenofile_and_get_samples.out.pheno_no_transform
        .map{["notransform", it]}
    }

    create_ancestry_x_transform_pheno(ch_transformed_phenos,
                                      ch_pca_keep_files)
    
    if (params.run_pca) {
      ch_input_for_add_pcs = filter_pca.out.pca_results.combine(create_ancestry_x_transform_pheno.out.crossed_pheno_out, by:0)
      add_pcs_to_pheno(ch_input_for_add_pcs)
      ch_final_phenos = add_pcs_to_pheno.out.final_phenos
    }

    if (!params.run_pca) {
      ch_final_phenos = create_ancestry_x_transform_pheno.out.crossed_pheno_out
    }

    align_pheno_with_test_variant_plink(ch_final_phenos,
                                        merge_plink.out.plink_merged)

    ch_aligned_test_vars_out = align_pheno_with_test_variant_plink.out.align_pheno_with_test_variant_plink_out
      .map{ it[0..1] + [it[2].text] + it[3..6] }
      .branch {
          binary: it[2] == "binary"
          quant: it[2] == "quantitative"
      }

    ch_aligned_test_vars_out_binary = ch_aligned_test_vars_out.binary

    ch_aligned_test_vars_out_quant = ch_aligned_test_vars_out.quant


    if (ch_aligned_test_vars_out_binary.count() != 0) {
      filter_binary_missingness(ch_aligned_test_vars_out_binary)
      ch_aligned_test_variants_plink = filter_binary_missingness.out.aligned_test_vars_binary_hwe_filtered
    }
    if (ch_aligned_test_vars_out_binary.count() == 0) {
      ch_aligned_test_variants_plink = Channel.empty()
    }

    if (ch_aligned_test_vars_out_quant.count() != 0) {
      ch_aligned_test_variants_plink.mix(ch_aligned_test_vars_out_quant)
    }

    //convert2bgen(ch_aligned_test_variants_plink)
    convert2bgen(ch_aligned_test_variants_plink)

    regenie_step1_fit_model(convert2bgen.out.merged_bgen)
    regenie_step2_association_testing(regenie_step1_fit_model.out.inputs_for_regenie_step2)
    regenie_association = regenie_step2_association_testing.out.regenie_step2_assoc

    //munge_regenie(regenie_step2_association_testing.out.regenie_out)

    summ_stats = munge_regenie.out.summ_stats

  emit:
    regenie_association
    summ_stats
    merged_plink
}


workflow{
  /*--------------------------------------------------
    Channel setup
  ---------------------------------------------------*/
  if (params.input_folder_location) {
    ch_user_input_vcf = Channel.fromPath("${params.input_folder_location}/**${params.file_pattern}*.{${params.file_suffix},${params.index_suffix}}")
        .map { it -> [ get_chromosome(file(it).simpleName.minus(".${params.index_suffix}").minus(".${params.file_suffix}")), "s3:/"+it] }
        .groupTuple(by:0)
        .map { chr, files_pair -> [ chr, files_pair[0], files_pair[1] ] }
        .map { chr, vcf, index -> [ file(vcf).simpleName, chr, file(vcf), file(index) ] }
        .take( params.number_of_files_to_process )
  }

  // KING tool for ancestry inference reference data
  ch_king_reference_data = Channel
    .fromFilePairs("${params.king_reference_data}",size:3, flat : true)
    .ifEmpty { exit 1, "KING reference data PLINK files not found: ${params.king_reference_data}.\nPlease specify a valid --king_reference_data value. e.g. refdata/king_ref*.{bed,bim,fam}" }

  //Pheno transform file
  ch_input_pheno_transform = Channel.fromPath("${params.pheno_transform}")


  // LocusZoom resource files
  //ch_hg38_gff3 = Channel.fromPath("${params.gff3_for_locuszoom}")

  projectDir = workflow.projectDir
  //ch_ancestry_inference_Rscript = Channel.fromPath("${projectDir}/bin/Ancestry_Inference.R", followLinks: false)
  //ch_DTable_Rscript = Channel.fromPath("${projectDir}/bin/DTable.R", followLinks: false)
  //ch_concat_chroms_Rscript = Channel.fromPath("${projectDir}/bin/concat_chroms.R", followLinks: false)
  //ch_convert_output_Rscript = Channel.fromPath("${projectDir}/bin/convert_output.R", followLinks: false)
  //ch_gwas_report_Rmd = Channel.fromPath("${projectDir}/bin/gwas_report.Rmd", followLinks: false)
  //ch_logo_png = Channel.fromPath("${projectDir}/bin/logo.png", followLinks: false)
  //ch_manhattan_Rscript = Channel.fromPath("${projectDir}/bin/manhattan.R", followLinks: false)
  //ch_pca_outliers_Rscript = Channel.fromPath("${projectDir}/bin/pca_outliers.R", followLinks: false)
  //ch_qqplot_Rscript = Channel.fromPath("${projectDir}/bin/qqplot.R", followLinks: false)
  //ch_sanitise_Rscriptch_sanitise_Rscript = Channel.fromPath("${projectDir}/bin/sanitise.R", followLinks: false)
  //ch_style_css = Channel.fromPath("${projectDir}/bin/style.css", followLinks: false)
  //ch_subset_gwascat_Rscript = Channel.fromPath("${projectDir}/bin/subset_gwascat.R", followLinks: false)
  //ch_transform_gwas_catalogue_Rscript = Channel.fromPath("${projectDir}/bin/transform_gwas_catalogue.R", followLinks: false)
  //ch_hail_gwas_script = Channel.fromPath("${projectDir}/bin/hail_gwas.py", followLinks: false)
  //ch_plotLocusZoom_script =  Channel.fromPath("${projectDir}/bin/plotLocusZoom.R", followLinks: false)
  //ch_locusZoomFunctions_script =  Channel.fromPath("${projectDir}/bin/locusZoomFunctions.R", followLinks: false)


  ch_user_input_vcf = Channel
    .fromPath(params.genotype_files_list)
    .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.genotype_files_list}" }
    .splitCsv(skip:1)
    .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
    .take( params.number_of_files_to_process )
  
  ch_high_ld_regions = Channel
    .fromPath(params.high_LD_long_range_regions)
    .ifEmpty { exit 1, "Cannot find file containing long-range LD regions for exclusion : ${params.high_LD_long_range_regions}" }

  ch_gwas_cat = Channel
    .fromPath(params.gwas_cat)
    .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }

  if (params.ld_scores) {
    ch_ld_scores = Channel
        .fromPath(params.ld_scores)
        .ifEmpty { exit 1, "Cannot find file containing LD scores : ${params.ld_scores}" }
  }

  ch_pheno = Channel
    .fromPath(params.pheno_data)
    .ifEmpty { exit 1, "Cannot find phenotype file : ${params.pheno_data}" }


  lifebitai_gwas_vcf_regenie(ch_user_input_vcf, 
                              ch_king_reference_data,
                              ch_input_pheno_transform,
                              ch_high_ld_regions,
                              ch_gwas_cat,
                              ch_ld_scores,
                              ch_pheno)
}