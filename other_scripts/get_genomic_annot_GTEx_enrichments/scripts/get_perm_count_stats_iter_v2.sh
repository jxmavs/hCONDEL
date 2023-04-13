type_iter=$1
label_ids_file=$2
type_perm_stats_file_list_file=$3
perm_type_results_dir=$4

type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_random
while read -a array
do
	label=${array[0]}
	label_file_id=${array[1]}
	rm -f "${perm_type_results_dir}/${type_iter}/${type_header_id}_${label_file_id}_count_summary_stats.txt"
	
    #check to see if there exists multiple count columns to loop over, if there are, they will be specified in ${labelFileNamesFile}.count_additional_meta, otherwise use the first field
    if [ -s "${label_ids_file}.count_additional_meta.${label_file_id}" ]
    then
        cp "${label_ids_file}.count_additional_meta.${label_file_id}" "${label_ids_file}.count_additional_meta.${label_file_id}.temp"
    else
        echo -e "${label}\t1\t${label_file_id}" > "${label_ids_file}.count_additional_meta.${label_file_id}.temp"
    fi
    count_meta_info_file="${label_ids_file}.count_additional_meta.${label_file_id}.temp"

	while read -a array2
	do
		perm_iter_file=${array2[0]}
		
		#get counts from all columns, put into ${all_counts}
		all_counts=""
		while IFS=$'\t' read -a array3
		do
			label_field_num="${array3[1]}"
			#get count
			count=$(awk -v label_field_num=$((${label_field_num}+1)) -v label=${label} '($label_field_num==label) {print}' ${perm_iter_file} | wc -l | awk '{print $1}')
			all_counts=${all_counts}"\t"${count}
		done<${count_meta_info_file}
		
		#remove last tab
		all_counts=$(echo -e ${all_counts} | sed 's/^\t//') 	
		echo -e "${all_counts}" >> ${perm_type_results_dir}/${type_iter}/${type_header_id}_${label_file_id}_count_summary_stats.txt

	done<${type_perm_stats_file_list_file}

done<${label_ids_file}

#count hCONDELs in any repeat class
if [ "${type_iter}" = "repeat_class" ]
then
	label="${type_iter}_all"
	label_file_id="${type_iter}_all"
	rm -f ${perm_type_results_dir}/${type_iter}/${type_header_id}_${label_file_id}_count_summary_stats.txt
	while read -a array
	do
		perm_iter_file=${array[0]}

		#get count
		count=$(awk -v label=${label} '$2!="NA" {print}' ${perm_iter_file} | wc -l | awk '{print $1}')

		echo -e "${count}" >> ${perm_type_results_dir}/${type_iter}/${type_header_id}_${label_file_id}_count_summary_stats.txt

	done<${type_perm_stats_file_list_file}
fi

### get stats of real data ###

type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
rm -f ${perm_type_results_dir}/${type_iter}/${type_header_id}_count_summary_stats.txt

while read -a array
do
	label=${array[0]}
	label_file_id=${array[1]}
	
	count_meta_info_file="${label_ids_file}.count_additional_meta.${label_file_id}.temp"
	
	#get counts from all columns, put into ${all_counts}
    all_counts=""
    while IFS=$'\t' read -a array3
    do
		label_field_num="${array3[1]}"
		#get count
		count=$(awk -v label_field_num=$((${label_field_num}+1)) -v label=${label} '($label_field_num==label) {print}' ${perm_type_results_dir}/${type_header_id}_info.txt | wc -l | awk '{print $1}')
		all_counts=${all_counts}"\t"${count}	
	done<${count_meta_info_file}
	
	#remove last tab
    all_counts=$(echo -e ${all_counts} | sed 's/^\t//') 
	echo -e "${label}\t${all_counts}" >> ${perm_type_results_dir}/${type_iter}/${type_header_id}_count_summary_stats.txt

	rm -f ${count_meta_info_file}
done<${label_ids_file}

#count hCONDELs in any repeat class
if [ "${type_iter}" = "repeat_class" ]
then
	count=$(awk -v label=${label} '($2!="NA") {print}' ${perm_type_results_dir}/${type_header_id}_info.txt | wc -l | awk '{print $1}')
	echo -e "${type_iter}_all\t${count}" >> ${perm_type_results_dir}/${type_iter}/${type_header_id}_count_summary_stats.txt
fi
