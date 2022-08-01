chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid
temp_dir=/cluster_path_temp/sgdp/temp
sgdp_output_dir=/broad/hptmp/sreilly/James
sgdp_output_dir=/cluster_path_temp/sgdp
af_filter=1
af_filter=0.00001
mkdir -p ${temp_dir}

rm -f ${temp_dir}/all_vcf_analysis.txt
while read -a array
do
    #this script analyzes the vcf from small deletions, and outputs preliminary statistics
    genome_file_id=${array[0]}
    echo "bash analyze_vcf_stats_post_call_new.sh ${genome_file_id} ${af_filter} ${sgdp_output_dir}">> ${temp_dir}/all_vcf_analysis.txt

done<${chimp_ref_file_names}




del_dir=/cluster_path/ape_project/deletions_project
totalMem=25G
totalNumTasks=$(wc -l ${temp_dir}/all_vcf_analysis.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_analyze_all_vcf.sh
qsub ${temp_dir}/submit_analyze_all_vcf.sh ${temp_dir}/all_vcf_analysis.txt
