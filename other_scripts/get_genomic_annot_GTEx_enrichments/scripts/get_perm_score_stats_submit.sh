type_ref_file_score=$1
perm_type_results_dir=$2
script_dir=$3

rm -f ${script_dir}/get_perm_score_stats_all.txt 
while read -a array
do
	type_iter=${array[0]}
	type_perm_stats_file_list_file=${array[3]}
	
	echo "bash '${script_dir}/get_perm_score_stats_iter.sh' '${type_iter}' '${type_perm_stats_file_list_file}' '${perm_type_results_dir}'" >> ${script_dir}/get_perm_score_stats_all.txt 
done<${type_ref_file_score}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=2G
totalNumTasks=$(wc -l ${script_dir}/get_perm_score_stats_all.txt  | awk '{print $1}')   
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/get_perm_score_stats_all.sh
totalTime=5:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/get_perm_score_stats_all.sh
rm -f ${script_dir}/get_perm_score_stats_all.log
qsub -o ${script_dir}/get_perm_score_stats_all.log ${script_dir}/get_perm_score_stats_all.sh ${script_dir}/get_perm_score_stats_all.txt

