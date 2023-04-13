barcode_master_table=$1
file_names_file=$2
seq_names_file=$3
temp_dir=$4
barcode_length=$5
header_name=$6
seq_id_filter_pct=$7
#############################################################################################################
script_path=/cluster_path/ape_project/deletions_project/barcode_associate/scripts

mkdir -p ${temp_dir}

rm -f ${temp_dir}/barcode_assoc_all_runs_${header_name}.txt
#submit jobs
while read -a array
do

    input_file_name="${array[0]}"
    output_file_name="${array[1]}"
    echo "bash ${script_path}/get_seq_table_counts_from_barcodes_v2.sh ${input_file_name} ${barcode_master_table} ${seq_names_file} ${barcode_length} ${seq_id_filter_pct} ${output_file_name}" >> ${temp_dir}/barcode_assoc_all_runs_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=150G
totalNumTasks=$(wc -l ${temp_dir}/barcode_assoc_all_runs_${header_name}.txt  | awk '{print $1}')
totalTime=100:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/barcode_assoc_submit_all_${header_name}.sh
rm -f ${temp_dir}/barcode_assoc_submit_all_${header_name}.log
qsub -o ${temp_dir}/barcode_assoc_submit_all_${header_name}.log  ${temp_dir}/barcode_assoc_submit_all_${header_name}.sh ${temp_dir}/barcode_assoc_all_runs_${header_name}.txt
