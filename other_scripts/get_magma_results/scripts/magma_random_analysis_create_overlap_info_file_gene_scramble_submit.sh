#!/bin/bash
#$ -cwd 
#$ -l h_vmem=3G
#$ -l h_rt=30:00:00
#$ -V
#$ -j y

random_analysis_header=$1
script_dir=$2
magma_random_analyses_type_temp_dir=$3
all_genes_vec_file_path=$4
hcondel_closest_dist_and_window_metadata_file_path=$5
closest_dist_filter=$6
window_size=$7
additional_info=$8 
output_file_path_with_header=$9
write_observed_overlap_genes=${10}
batch_info_file=${11}
write_random_overlap_genes=${12}

rm -f ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.txt

while IFS=$'\t' read -a array
do
	batch_num="${array[0]}"
	
	additional_info_iter=${additional_info}_batch_${batch_num}

	#remove any old files
    rm -f ${output_file_path_with_header}_closest_dist_${closest_dist_filter}_${additional_info_iter}_hg38_protein_coding_gene_overlap_info.txt
    rm -f ${output_file_path_with_header}_closest_dist_${closest_dist_filter}_${additional_info_iter}_random_matched_hg38_protein_coding_gene_overlap_info.txt 
	
	rm -f ${output_file_path_with_header}_window_${window_size}_${additional_info_iter}_hg38_protein_coding_gene_overlap_info.txt
	rm -f ${output_file_path_with_header}_window_${window_size}_${additional_info_iter}_random_matched_hg38_protein_coding_gene_overlap_info.txt
	
	#generate set_file from magma_analysis_preprocess_2_v2.r 
	#first column gene id, second column group
	echo -e "Rscript ${script_dir}/magma_analysis_preprocess_2_v2.r ${all_genes_vec_file_path} ${hcondel_closest_dist_and_window_metadata_file_path} ${closest_dist_filter} ${window_size} ${additional_info_iter} ${output_file_path_with_header} ${write_observed_overlap_genes} ${batch_num} ${write_random_overlap_genes}" >> ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.txt
 	
done<${batch_info_file}

#qsub parameters
totalMem=3G
totalNumTasks=$(wc -l ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.txt | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${script_dir}/runTasksEval.sh > ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.sh

rm -f ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.log
qsub -o ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.log ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.sh ${magma_random_analyses_type_temp_dir}/${random_analysis_header}_create_overlap_info_file.txt

