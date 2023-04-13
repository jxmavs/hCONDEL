sgdp_file_name=$1
sgdp_output_dir=$2
chimp_ref_file_names=$3

#align sgdp unitigs
sgdp_dir=/cluster_path_temp/sgdp
#ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_seq_1
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_2.txt

chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid
sgdp_indiv_names_file=${sgdp_dir}/${sgdp_file_name}
temp_dir=/cluster_path_temp/sgdp/temp

rm -f ${temp_dir}/all_runs_call_variants.txt
 
while read -a array
do
    genome_file_id=${array[0]}
   
    sgdp_bam_dir=${sgdp_output_dir}/${genome_file_id}_bam_files
    sgdp_vcf_dir=${sgdp_output_dir}/${genome_file_id}_vcf_files
    mkdir -p ${sgdp_bam_dir}
    mkdir -p ${sgdp_vcf_dir}
    
    #sgdp_sam_dir=/cluster_path_temp/sgdp/${genome_file_id}_sam_files
    #ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_seq_1.fa
    chimp_human_hybrid_out_fasta_name=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa
    echo -e "bash sgdp_call_chimp_hybrid_variants_iter.sh ${sgdp_bam_dir} ${sgdp_vcf_dir} ${chimp_human_hybrid_out_fasta_name} ${sgdp_indiv_names_file} ${genome_file_id}" >> ${temp_dir}/all_runs_call_variants.txt

done<${chimp_ref_file_names}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=3G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_call_variants.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_all_call_variants.sh
qsub ${temp_dir}/submit_all_call_variants.sh ${temp_dir}/all_runs_call_variants.txt

