#align sgdp unitigs

sgdp_dir=/cluster_path_temp/sgdp
#ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_seq_1
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_2.txt
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid
del_dir=/cluster_path/ape_project/deletions_project 
temp_dir=/cluster_path_temp/sgdp/temp
#mkdir -p ${temp_dir}
#temp_dir=$1
sgdp_input_dir=$1
sgdp_output_dir=$2
sgdp_file_name=$3
chimp_ref_file_names=$4
numcores=1
 
while read -a array
do
    #echo ${sgdp_dir}/${sgdp_file_name}
    genome_file_id=${array[0]}
    ref_file_path=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord
    rm -f ${temp_dir}/all_runs_align_${genome_file_id}.txt
    #put files in bam directory
    sgdp_output_genome_file_id_dir=${sgdp_output_dir}/${genome_file_id}_bam_files
    mkdir -p ${sgdp_output_genome_file_id_dir}
    while read -a array
    do
    
        pop_sample_file_path=${sgdp_input_dir}/${array[0]}
        pop_sample_id=${array[1]}
        echo -e "bash ${del_dir}/bwa_align_tigs_iter.sh ${sgdp_input_dir} ${sgdp_output_genome_file_id_dir} ${ref_file_path} ${pop_sample_file_path} ${pop_sample_id} ${genome_file_id} ${numcores}" >> ${temp_dir}/all_runs_align_${genome_file_id}.txt
        #echo ${pop_sample_file_path}
        #echo ${pop_sample_id}
    done<${sgdp_dir}/${sgdp_file_name}

    #submit all jobs
    totalMem=15G
    #echo ${temp_dir}/all_runs_align_${genome_file_id}.txt
    totalNumTasks=$( wc -l ${temp_dir}/all_runs_align_${genome_file_id}.txt  | awk '{print $1}' )
    sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -q short/#$ -q long/g; s/#$ -pe smp numCores/#$ -pe smp ${numcores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numcores}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/submit_all_align_${genome_file_id}.sh
    #sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -pe smp numCores/#$ -pe smp ${numcores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numcores}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/submit_all_align_${genome_file_id}.sh
    qsub -o ${temp_dir}/submit_all_align_${genome_file_id}.log ${temp_dir}/submit_all_align_${genome_file_id}.sh ${temp_dir}/all_runs_align_${genome_file_id}.txt

done<${chimp_ref_file_names}
