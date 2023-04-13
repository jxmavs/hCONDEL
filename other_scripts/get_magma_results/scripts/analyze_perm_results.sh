#get summary stats
analysis_run_name=$1
column_metadata_file=$2
empirical_stats_file=$3
perm_stats_file=$4
perm_analysis_out_file=$5

rm -f ${perm_analysis_out_file}

while IFS=$'\t' read -a array
do 
	column_number="${array[0]}"
	column_label="${array[1]}"
	test_direction="${array[2]}"
		
	empirical_value=$(awk -F"\t" -v column_number=${column_number} '{print $column_number}' ${empirical_stats_file})
	
	if [ "${test_direction}" = "lt" ]
	then
		numerator_value=$(awk -F"\t" -v column_number=${column_number} -v empirical_value=${empirical_value} '$column_number<empirical_value{print}' ${perm_stats_file} | wc -l | awk '{print $1}')
	else
		numerator_value=$(awk -F"\t" -v column_number=${column_number} -v empirical_value=${empirical_value} '$column_number>empirical_value{print}' ${perm_stats_file} | wc -l | awk '{print $1}')	
	fi
	
	denominator_value=$(awk '{print}' ${perm_stats_file} | wc -l | awk '{print $1}')
	
	if [ ${denominator_value} -eq 0 ]
	then 
		p_value="NA"
	else
		p_value=$(echo -e "${numerator_value}\t${denominator_value}" | awk -F"\t" '{print $1/$2}')    
	fi

	perm_column_mean=$(awk -v column_number=${column_number} '{ sum += $column_number; n++ } END { if (n > 0) print sum / n; }' ${perm_stats_file})
	perm_column_sd=$(awk -v column_number=${column_number} '{sum += $column_number; sumsq += ($column_number)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_stats_file})

	#https://stackoverflow.com/questions/2424770/floating-point-comparison-in-shell
	var=$(awk 'BEGIN{ print "'${perm_column_sd}'"==0 }') 
	if [ ${var} -eq 1 ]
	then
		z_score="NA"
	else
		z_score=$( echo -e "${empirical_value}\t${perm_column_mean}\t${perm_column_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
	fi
	
	echo -e "${analysis_run_name}\t${column_label}\t${empirical_value}\t${p_value}\t${perm_column_mean}\t${perm_column_sd}\t${z_score}" >> ${perm_analysis_out_file}
done<${column_metadata_file}


