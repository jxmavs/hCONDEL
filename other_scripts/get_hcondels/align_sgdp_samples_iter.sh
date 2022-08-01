
#align sgdp unitigs
sgdp_dir=/cluster_path_temp/sgdp
sgdp_file_name=sgdp_file_names_2.txt
#ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_seq_1

ref_file_path=$1
identifier=$2
while read -a array
do
    
    pop_sample_file_name=${array[0]}
    pop_sample_id=${array[1]}
    echo -e "bash bwa_align_tigs_iter.sh ${sgdp_dir} ${ref_file_path} ${pop_sample_file_name} ${pop_sample_id} ${identifier}"
    #echo ${pop_sample_file_name}
    #echo ${pop_sample_id}
done<${sgdp_dir}/${sgdp_file_name}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=3G
totalNumTasks=$(wc -l ${temp_dir}/all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_all.sh
qsub ${temp_dir}/submit_all.sh ${temp_dir}/all_runs.txt
