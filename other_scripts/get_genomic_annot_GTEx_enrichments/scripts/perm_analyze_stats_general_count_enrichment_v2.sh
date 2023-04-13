#get summary stats
type_iter=$1
label_ids_file=$2
real_stats_file=$3
outFile=$4
perm_summary_stats_file_header=$5

rm -f ${outFile}

if [ "${type_iter}" = "repeat_class" ]
then
    cat ${label_ids_file} > ${label_ids_file}.temp
    echo -e "${type_iter}_all\t${type_iter}_all" >> ${label_ids_file}.temp
else
    cat ${label_ids_file} > ${label_ids_file}.temp
fi
#cat ${label_ids_file} > ${label_ids_file}.temp

while IFS=$'\t' read -a array
do
	label=${array[0]}
	label_file_id=${array[1]}
	perm_summary_stats_file="${perm_summary_stats_file_header}_${label_file_id}_count_summary_stats.txt"    
	
	#check to see if there exists multiple count columns to loop over, if there are, they will be specified in ${label_ids_file}.count_additional_meta, otherwise use the first field
	if [ -s "${label_ids_file}.count_additional_meta.${label_file_id}" ]
	then
		cp "${label_ids_file}.count_additional_meta.${label_file_id}" "${label_ids_file}.count_additional_meta.${label_file_id}.temp"
	else
		echo -e "${label}\t1\t${type_iter}" > "${label_ids_file}.count_additional_meta.${label_file_id}.temp"
	fi
	count_meta_info_file="${label_ids_file}.count_additional_meta.${label_file_id}.temp"
	
	while IFS=$'\t' read -a array2
	do 
		label="${array2[0]}"
		label_field_num="${array2[1]}"
		label_out="${array2[2]}"
		
		######### for count ##############
		#get real value
		#first field is label in real stats file, so add one to correct
		real_count=$(awk -F"\t" -v label_field_num=$((${label_field_num}+1)) -v label=${label} '$1==label {print $label_field_num}' ${real_stats_file})
		
		count_num=$(awk -F"\t" -v label_field_num=${label_field_num} -v real_count=${real_count} '$label_field_num>real_count{print}' ${perm_summary_stats_file} | wc -l)
		count_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)
		
		if [ ${count_denom} -eq 0 ]
		then 
			count_pval="NA"
		else
			count_pval=$(echo -e "${count_num}\t${count_denom}" | awk -F"\t" '{print $1/$2}')    
		fi

		perm_count_mean=$(awk -v label_field_num=${label_field_num} '{ sum += $label_field_num; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
		perm_count_sd=$(awk -v label_field_num=${label_field_num} '{sum += $label_field_num; sumsq += ($label_field_num)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})

		#https://stackoverflow.com/questions/2424770/floating-point-comparison-in-shell
		var=$(awk 'BEGIN{ print "'${perm_count_sd}'"==0 }') 
		if [ ${var} -eq 1 ]
		then
			count_z_score="NA"
		else
			count_z_score=$( echo -e "${real_count}\t${perm_count_mean}\t${perm_count_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
		fi
		
		#echo -e "${label}\t${real_count}\t${count_pval}\t${perm_count_mean}\t${perm_count_sd}\t${count_z_score}"
		echo -e "${label}\t${real_count}\t${count_pval}\t${perm_count_mean}\t${perm_count_sd}\t${count_z_score}\t${label_out}" >> ${outFile}
	done<${count_meta_info_file}
	
	rm -f ${count_meta_info_file}
done<${label_ids_file}.temp

rm -f ${label_ids_file}.temp
