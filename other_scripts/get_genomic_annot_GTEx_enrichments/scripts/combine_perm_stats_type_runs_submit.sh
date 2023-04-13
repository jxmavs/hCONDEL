type_ref_file=$1
perm_type_results_dir=$2
script_dir=$3

rm -f ${script_dir}/combine_perm_stats_type_runs.txt 
while read -a array
do
	type_iter=${array[0]}
	
	echo "bash '${script_dir}/combine_perm_stats_type_runs_iter.sh' '${type_iter}' '${perm_type_results_dir}'" >> ${script_dir}/combine_perm_stats_type_runs.txt 
done<${type_ref_file}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=2G
totalNumTasks=$(wc -l ${script_dir}/combine_perm_stats_type_runs.txt  | awk '{print $1}')   
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/combine_perm_stats_type_runs.sh
totalTime=3:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${script_dir}/combine_perm_stats_type_runs.sh
rm -f ${script_dir}/combine_perm_stats_type_runs.log
qsub -o ${script_dir}/combine_perm_stats_type_runs.log ${script_dir}/combine_perm_stats_type_runs.sh ${script_dir}/combine_perm_stats_type_runs.txt

