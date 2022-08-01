
#this script creates the chimp human reference genomes for each coordinate
human_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human
filter_dir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

temp_dir=/cluster_path_temp/sgdp/temp
mkdir -p ${temp_dir}
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
mkdir -p ${temp_dir}


rm -f ${temp_dir}/all_runs_build_hybrid_index.txt
#submit jobs
while read -a array
do
    genome_file_id=${array[0]}
    #chimp_human_hybrid_out_fasta_name="${array[3]}"
    #chimp_human_hybrid_out_index_file_name=$(echo ${chimp_human_hybrid_out_fasta_name} | awk -F'.' '{print $(NF-1)}') 
    chimp_human_hybrid_out_fasta_name=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa
    chimp_human_hybrid_out_index_file_name=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord 
    echo "bash build_chimp_hybrid_index_iter.sh ${chimp_human_hybrid_out_index_file_name} ${chimp_human_hybrid_out_fasta_name}" >> ${temp_dir}/all_runs_build_hybrid_index.txt
    
done<${chimp_ref_file_names}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project 
totalMem=10G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_build_hybrid_index.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_build_hybrid_index.sh
qsub ${temp_dir}/submit_build_hybrid_index.sh ${temp_dir}/all_runs_build_hybrid_index.txt

