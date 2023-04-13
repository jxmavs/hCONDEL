#combine all files from fimo run into one master file
out_header=$1
#the file header for the resulting output
file_out_header=$2
out_path=$3
pairwise_comp_file=$4
tf_name=$5
seqFileNamesFile=$6

############################

while read -a array
do
	
	input_file=${array[0]}
	file_id_suffix=$( echo ${input_file} | awk -F"." '{print $NF}' )
	suffix=${file_id_suffix}"."${tf_name}
	
	rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output.${suffix}.txt

	if [[ "${out_header}" == "tf_enrichment_hcondel_denylist_filtered_random_JASPARv2020" ]]
	then 
		#rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.${tf_name}.txt	
		rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt
		rm -fr ${out_path}/perm_final/${tf_name}
		rm -f ${out_path}/perm_final/${tf_name}_file_list.txt
		rm -f ${out_path}/perm_final/${tf_name}_file_check.txt
	fi
	
done<${seqFileNamesFile}



