out_path=$1
temp_dir=$2
script_dir=$3
header_name=$4
seq_table_file_path=$5
num_jobs_per_task=$6

#make sure all files exist

fileInd=0
taskFileInd=0

run_name=get_mismatch_stats_v2_mult_jobs_per_task

rm -f ${temp_dir}/${run_name}_${header_name}.txt


split_tasks_submit_run_file=${out_path}/get_mismatch_stats_v2_split_tasks_submit.sh
rm -f ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh
rm -f ${split_tasks_submit_run_file}

while read -a array
do
    seq_name="${array[0]}"
    chimp_seq="${array[1]}"
    human_seq="${array[2]}"
    seq_id="${array[3]}"
     
    echo -e "bash ${script_dir}/get_mismatch_stats_iter_v2.sh ${out_path} ${header_name} \"${seq_name}\" ${chimp_seq} ${human_seq} ${seq_id}" >> ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh

    fileInd=$((${fileInd}+1))

    #if we have hit the num_jobs_per_task, then move to the next task file
    if [[ $(( ${fileInd}%${num_jobs_per_task} )) == 0 ]]; then
        echo -e "bash ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh" >> ${split_tasks_submit_run_file}
        taskFileInd=$((${taskFileInd}+1))
        rm -f ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh
    fi
    
done<${seq_table_file_path}

#put in any remaining tasks
if [[ $(( ${fileInd}%${num_jobs_per_task} )) != 0 ]]; then
    echo -e "bash ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh" >> ${split_tasks_submit_run_file}
    taskFileInd=$((${taskFileInd}+1))
    rm -f ${out_path}/${run_name}_split_tasks_submit_task_${taskFileInd}.sh
fi




#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=1G
totalNumTasks=$(wc -l ${split_tasks_submit_run_file} | awk '{print $1}')
totalTime=50:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/${run_name}_${header_name}.sh
rm -f ${temp_dir}/${run_name}_${header_name}.log
qsub -o ${temp_dir}/${run_name}_${header_name}.log ${temp_dir}/${run_name}_${header_name}.sh ${split_tasks_submit_run_file} 
