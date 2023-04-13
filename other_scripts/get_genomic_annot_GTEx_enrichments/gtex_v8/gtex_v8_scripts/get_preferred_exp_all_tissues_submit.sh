sample_gene_exp_table_file=$1
sample_attributes_table_file=$2
out_file_path_with_file_header=$3
tissue_list_file=$4

script_dir=/cluster_path/ape_project/deletions_project/gtex_v8/scripts

rm -f ${out_file_path_with_file_header}_get_preferred_exp_all_runs.txt 
while IFS=$'\t' read -a array
do
	tissue_of_interest="${array[0]}"
	echo "Rscript ${script_dir}/get_tissue_specific_genes.R ${sample_gene_exp_table_file} ${sample_attributes_table_file} ${out_file_path_with_file_header} '${tissue_of_interest}'" >> ${out_file_path_with_file_header}_get_preferred_exp_all_runs.txt 
done<${tissue_list_file}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=60G
totalNumTasks=$(wc -l ${out_file_path_with_file_header}_get_preferred_exp_all_runs.txt  | awk '{print $1}')   
totalTime=50:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${out_file_path_with_file_header}_get_preferred_exp.sh
rm -f ${out_file_path_with_file_header}_get_preferred_exp.log
qsub -o ${out_file_path_with_file_header}_get_preferred_exp.log ${out_file_path_with_file_header}_get_preferred_exp.sh ${out_file_path_with_file_header}_get_preferred_exp_all_runs.txt

