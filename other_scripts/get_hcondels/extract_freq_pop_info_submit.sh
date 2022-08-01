af_filter=0.00001
sgdp_output_dir=/cluster_path_temp/sgdp 
del_dir=/cluster_path/ape_project/deletions_project 
#sgdp_meta_file is extracted from the sgdp paper
sgdp_meta_file=${del_dir}/sgdp_region_country_id_cleaned.txt
temp_dir=/cluster_path_temp/sgdp/temp
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt

#this file extracts the polymorphic indels and obtains info on the allele frequency, the populations that they are varying, etc.
#it should be run after analyze_vcf_stats_submit.sh
rm -f ${temp_dir}/all_vcf_freq_analysis.txt
while read -a array 
do
    genome_file_id=${array[0]}
    vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
    echo "bash extract_freq_pop_info_iter.sh ${genome_file_id} ${af_filter} ${sgdp_output_dir} ${sgdp_meta_file}">> ${temp_dir}/all_vcf_freq_analysis.txt   

done<${chimp_ref_file_names}

del_dir=/cluster_path/ape_project/deletions_project
totalMem=10G
totalNumTasks=$(wc -l ${temp_dir}/all_vcf_freq_analysis.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_analyze_freq_all_vcf.sh
qsub ${temp_dir}/submit_analyze_freq_all_vcf.sh ${temp_dir}/all_vcf_freq_analysis.txt
