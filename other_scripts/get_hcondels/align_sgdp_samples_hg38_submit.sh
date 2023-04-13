#align sgdp unitigs

sgdp_dir=/cluster_path_temp/sgdp
#ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_seq_1
del_dir=/cluster_path/ape_project/deletions_project 
temp_dir=/cluster_path_temp/sgdp/temp
#mkdir -p ${temp_dir}
sgdp_input_dir=$1
sgdp_output_dir=$2
sgdp_file_name=$3
ref_file_path=$4
genome_file_id=$5
numcores=1

#echo ${sgdp_dir}/${sgdp_file_name}
rm -f ${temp_dir}/all_runs_align_${genome_file_id}.txt
#put files in bam directory
sgdp_output_genome_file_id_dir=${sgdp_output_dir}/${genome_file_id}_bam_files
mkdir -p ${sgdp_output_genome_file_id_dir}
while read -a array
do

    pop_sample_file_path=${sgdp_input_dir}/${array[0]}
    pop_sample_id=${array[1]}
    echo -e "bash ${del_dir}/bwa_align_tigs_iter.sh ${sgdp_output_genome_file_id_dir} ${ref_file_path} ${pop_sample_file_path} ${pop_sample_id} ${genome_file_id} ${numcores}" >> ${temp_dir}/all_runs_align_${genome_file_id}.txt
    #echo ${pop_sample_file_path}
    #echo ${pop_sample_id}
done<${sgdp_dir}/${sgdp_file_name}

#submit all jobs
totalMem=15G
#echo ${temp_dir}/all_runs_align_${genome_file_id}.txt
totalNumTasks=$( wc -l ${temp_dir}/all_runs_align_${genome_file_id}.txt  | awk '{print $1}' )
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -q short/#$ -q long/g; s/#$ -pe smp numCores/#$ -pe smp ${numcores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numcores}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/submit_all_align_${genome_file_id}.sh
#sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -pe smp numCores/#$ -pe smp ${numcores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numcores}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/submit_all_align_${genome_file_id}.sh
rm -f ${temp_dir}/submit_all_align_${genome_file_id}.log 
qsub -o ${temp_dir}/submit_all_align_${genome_file_id}.log ${temp_dir}/submit_all_align_${genome_file_id}.sh ${temp_dir}/all_runs_align_${genome_file_id}.txt

