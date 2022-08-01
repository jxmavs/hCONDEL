
#this script creates the chimp human reference genomes for each coordinate
human_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human
filter_dir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

temp_dir=/cluster_path_temp/sgdp/temp
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
mkdir -p ${temp_dir}

human_fasta_file_name=${human_dir}/hg38.fa 
human_fasta_chr_sizes=${human_dir}/hg38.chrom.sizes

rm -f ${temp_dir}/all_runs_create_hybrid_genome.txt
#submit jobs
while read -a array
do
    file_name=${array[0]}
    #chimp_coord_file_name="${array[0]}"
    #chimp_coord_fasta_name="${array[1]}"
    #chimp_human_hybrid_out_file_name="${array[2]}"
    #chimp_human_hybrid_out_meta_name="${array[3]}"
    chimp_coord_file_name=${chimp_human_hybrid_dir}/${file_name}.txt
    chimp_coord_fasta_name=${chimp_human_hybrid_dir}/${file_name}.fa
    chimp_human_hybrid_out_meta_name=${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.txt
    chimp_human_hybrid_out_fasta_name=${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.fa
    echo "bash create_hybrid_sequence_iter.sh ${chimp_coord_file_name} ${chimp_coord_fasta_name} ${human_fasta_file_name} ${human_fasta_chr_sizes} ${chimp_human_hybrid_out_fasta_name} ${chimp_human_hybrid_out_meta_name}" >> ${temp_dir}/all_runs_create_hybrid_genome.txt
    
done<${chimp_ref_file_names}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=30G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_create_hybrid_genome.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_all_create_hybrid_genome.sh
qsub ${temp_dir}/submit_all_create_hybrid_genome.sh ${temp_dir}/all_runs_create_hybrid_genome.txt

