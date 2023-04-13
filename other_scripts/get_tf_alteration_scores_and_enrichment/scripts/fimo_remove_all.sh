
out_path=$1

cell_type_file=/cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_cell_list.txt
tf_db_file=/cluster_path/ape_project/deletions_project/tf_misc_files/tf_db_file.txt

while read -a array
do
	cell_header_id=${array[0]}
	while read -a array1
	do
		tf_db_name=${array1[0]}
		out_header=${cell_header_id}_del_encode_annotation_${tf_db_name}
		while read -a array2
		do
			
			input_file=${array2[0]}
			suffix=$( echo ${input_file} | awk -F"." '{print $NF}' )
			rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output.${suffix}.txt 

			rm -f ${out_path}/${out_header}_fimo_cleaned_sum_diff_output.${suffix}.txt 

		done<${out_path}/all_fimo_file_names.txt
	done<${tf_db_file}
done<${cell_type_file}

