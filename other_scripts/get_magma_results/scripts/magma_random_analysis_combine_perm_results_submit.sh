#!/bin/bash
#$ -cwd 
#$ -l h_vmem=3G
#$ -l h_rt=30:00:00
#$ -V
#$ -j y

random_analysis_header=$1
gwas_analyses_of_interest_metadata_file=$2
batch_info_file=$3
script_dir=$4
magma_random_analyses_type_dir=$5
magma_random_analyses_type_temp_dir=$6

rm -f ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.txt

while IFS=$'\t' read -a array
do

    gwas_analysis_header="${array[0]}"

	#remove old files
	rm -f ${magma_random_analyses_type_dir}/${gwas_analysis_header}/${random_analysis_header}_${gwas_analysis_header}_all_perm_stats.txt
		
	echo -e "bash ${script_dir}/magma_random_analysis_combine_perm_results.sh ${random_analysis_header} ${gwas_analysis_header} ${gwas_analyses_of_interest_metadata_file} ${batch_info_file} ${magma_random_analyses_type_dir}" >> ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.txt

done<${gwas_analyses_of_interest_metadata_file}

#qsub parameters
totalMem=3G
totalNumTasks=$(wc -l ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.txt | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${script_dir}/runTasksEval.sh > ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.sh

rm -f ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.log
qsub -o ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.log ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.sh ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_combine_perm_results.txt

