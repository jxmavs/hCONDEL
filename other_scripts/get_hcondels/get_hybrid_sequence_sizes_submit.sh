
#this script gets the sizes for the hybrid genome
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

temp_dir=/cluster_path_temp/sgdp/temp
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
mkdir -p ${temp_dir}


rm -f ${temp_dir}/all_runs_create_hybrid_genome.txt
#submit jobs
while read -a array
do
    file_name=${array[0]}
    chimp_human_hybrid_fasta_name=${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.fa
    chimp_human_hybrid_fasta_sizes=${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.sizes
    echo "bash get_hybrid_sequence_sizes_iter.sh ${chimp_human_hybrid_fasta_name} ${chimp_human_hybrid_fasta_sizes}" >> ${temp_dir}/all_runs_get_hybrid_genome_sizes.txt
    
done<${chimp_ref_file_names}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=30G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_get_hybrid_genome_sizes.txt | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_all_get_hybrid_genome_sizes.sh 
qsub ${temp_dir}/submit_all_get_hybrid_genome_sizes.sh ${temp_dir}/all_runs_get_hybrid_genome_sizes.txt

