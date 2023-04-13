input_file_name=$1
barcode_master_table=$2
seq_names_file=$3
barcode_length=$4
seq_id_filter_pct=$5
output_file_name=$6

script_path=/cluster_path/ape_project/deletions_project/barcode_associate/scripts
python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin

#gunzip -c ${input_file_name}.gz > ${input_file_name}
#take gzip file straight
${python_dir}/python ${script_path}/barcodeToDelSeqV2.py ${input_file_name}.gz ${barcode_master_table} ${seq_names_file} ${barcode_length} ${seq_id_filter_pct} ${output_file_name}

#rm -f ${input_file_name}
