type_ref_file=$1
perm_type_results_dir=$2
script_dir=$3

rm -f ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.txt 
while read -a array
do
    type_iter=${array[0]}
    label_ids_file=${array[2]}

    type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
    real_stats_file=${perm_type_results_dir}/${type_iter}/${type_header_id}_count_summary_stats.txt

    type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_random
    perm_summary_stats_file_header=${perm_type_results_dir}/${type_iter}/${type_header_id}

    out_path=${perm_type_results_dir}/${type_iter}
    outFile=${out_path}/${type_iter}_count_enrichment_hcondel_denylist_filtered_final_stats.txt

    echo -e "bash ${script_dir}/perm_analyze_stats_general_count_enrichment_v2.sh '${type_iter}' '${label_ids_file}' '${real_stats_file}' '${outFile}' '${perm_summary_stats_file_header}'" >> ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.txt
done<${type_ref_file}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=2G
totalNumTasks=$(wc -l ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.txt  | awk '{print $1}')   
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.sh
totalTime=3:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.sh
rm -f ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.log
qsub -o ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.log ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.sh ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.txt

