#get summary stats

memeFileNamesFile=$1
real_stats_file=$2
outFile=$3
perm_summary_stats_file_header=$4

rm -f ${outFile}

echo -e "tf_name\tempirical_count\tcount_pval\tperm_count_mean\tperm_count_sd\tcount_z_score\tempirical_binding_diff\tbinding_diff_pval\tperm_binding_diff_mean\tperm_binding_diff_sd\tbinding_diff_z_score\tempirical_count_strengthen\tcount_strengthen_pval\tperm_count_strengthen_mean\tperm_count_strengthen_sd\tcount_strengthen_z_score\tempirical_binding_diff_strengthen\tbinding_diff_strengthen_pval\tperm_binding_diff_strengthen_mean\tperm_binding_diff_strengthen_sd\tbinding_diff_strengthen_z_score\tempirical_count_weaken\tcount_weaken_pval\tperm_count_weaken_mean\tperm_count_weaken_sd\tcount_weaken_z_score\tempirical_binding_diff_weaken\tbinding_diff_weaken_pval\tperm_binding_diff_weaken_mean\tperm_binding_diff_weaken_sd\tbinding_diff_weaken_z_score" > ${outFile}

while IFS=$'\t' read -a array
do
	motif_file=${array[0]}
    tf_name=$( echo ${motif_file} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' )

	perm_summary_stats_file=${perm_summary_stats_file_header}_${tf_name}_summary_stats.txt

    ######### for count ##############
    #get real value
    real_count=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $2}' ${real_stats_file})
    
	count_num=$(awk -F"\t" -v real_count=${real_count} '$1>real_count{print}' ${perm_summary_stats_file} | wc -l)
    count_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)
	count_pval=$(echo -e "${count_num}\t${count_denom}" | awk -F"\t" '{print $1/$2}')

    perm_count_mean=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
	perm_count_sd=$(awk '{sum += $1; sumsq += ($1)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})
	
	var=$(awk 'BEGIN{ print "'${perm_count_sd}'"==0 }')
	if [ ${var} -eq 1 ]
	then
		count_z_score="NA"
	else
		count_z_score=$( echo -e "${real_count}\t${perm_count_mean})\t${perm_count_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
	fi
	
	######### for binding_diff_abs_avg ##############
    #get real value
    real_binding_diff_abs_avg=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $3}' ${real_stats_file})
    
    binding_diff_abs_avg_num=$(awk -F"\t" -v real_binding_diff_abs_avg=${real_binding_diff_abs_avg} '$2>real_binding_diff_abs_avg{print}' ${perm_summary_stats_file} | wc -l)
    binding_diff_abs_avg_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)	
	binding_diff_abs_avg_pval=$(echo -e "${binding_diff_abs_avg_num}\t${binding_diff_abs_avg_denom}" | awk -F"\t" '{print $1/$2}')    
    
    perm_binding_diff_abs_avg_mean=$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_binding_diff_abs_avg_sd=$(awk '{sum += $2; sumsq += ($2)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})
    
    var=$(awk 'BEGIN{ print "'${perm_binding_diff_abs_avg_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
		binding_diff_abs_avg_z_score="NA"
	else
        binding_diff_abs_avg_z_score=$( echo -e "${real_binding_diff_abs_avg}\t${perm_binding_diff_abs_avg_mean})\t${perm_binding_diff_abs_avg_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi

    ######### for count_strengthen ##############
    #get real value
    real_count_strengthen=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $4}' ${real_stats_file})

    count_strengthen_num=$(awk -F"\t" -v real_count_strengthen=${real_count_strengthen} '$3>real_count_strengthen{print}' ${perm_summary_stats_file} | wc -l)
    count_strengthen_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)
    count_strengthen_pval=$(echo -e "${count_strengthen_num}\t${count_strengthen_denom}" | awk -F"\t" '{print $1/$2}')

    perm_count_strengthen_mean=$(awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_count_strengthen_sd=$(awk '{sum += $3; sumsq += ($3)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})

    var=$(awk 'BEGIN{ print "'${perm_count_strengthen_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
        count_strengthen_z_score="NA"
    else
        count_strengthen_z_score=$( echo -e "${real_count_strengthen}\t${perm_count_strengthen_mean})\t${perm_count_strengthen_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi

    ######### for binding_diff_strengthen_abs_avg ##############
    #get real value
    real_binding_diff_strengthen_abs_avg=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $5}' ${real_stats_file})
    
    binding_diff_strengthen_abs_avg_num=$(awk -F"\t" -v real_binding_diff_strengthen_abs_avg=${real_binding_diff_strengthen_abs_avg} '$4>real_binding_diff_strengthen_abs_avg{print}' ${perm_summary_stats_file} | wc -l)
    binding_diff_strengthen_abs_avg_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)  
    binding_diff_strengthen_abs_avg_pval=$(echo -e "${binding_diff_strengthen_abs_avg_num}\t${binding_diff_strengthen_abs_avg_denom}" | awk -F"\t" '{print $1/$2}')    
    
    perm_binding_diff_strengthen_abs_avg_mean=$(awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_binding_diff_strengthen_abs_avg_sd=$(awk '{sum += $4; sumsq += ($4)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})
    
    var=$(awk 'BEGIN{ print "'${perm_binding_diff_strengthen_abs_avg_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
        binding_diff_strengthen_abs_avg_z_score="NA"
    else
        binding_diff_strengthen_abs_avg_z_score=$( echo -e "${real_binding_diff_strengthen_abs_avg}\t${perm_binding_diff_strengthen_abs_avg_mean})\t${perm_binding_diff_strengthen_abs_avg_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi

    ######### for count_weaken ##############
    #get real value
    real_count_weaken=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $6}' ${real_stats_file})

    count_weaken_num=$(awk -F"\t" -v real_count_weaken=${real_count_weaken} '$5>real_count_weaken{print}' ${perm_summary_stats_file} | wc -l)
    count_weaken_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)
    count_weaken_pval=$(echo -e "${count_weaken_num}\t${count_weaken_denom}" | awk -F"\t" '{print $1/$2}')

    perm_count_weaken_mean=$(awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_count_weaken_sd=$(awk '{sum += $5; sumsq += ($5)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})

    var=$(awk 'BEGIN{ print "'${perm_count_weaken_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
        count_weaken_z_score="NA"
    else
        count_weaken_z_score=$( echo -e "${real_count_weaken}\t${perm_count_weaken_mean})\t${perm_count_weaken_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi

    ######### for binding_diff_weaken_abs_avg ##############
    #get real value
    real_binding_diff_weaken_abs_avg=$(awk -F"\t" -v tf_name=${tf_name} '$1==tf_name {print $7}' ${real_stats_file})
    
    binding_diff_weaken_abs_avg_num=$(awk -F"\t" -v real_binding_diff_weaken_abs_avg=${real_binding_diff_weaken_abs_avg} '$6>real_binding_diff_weaken_abs_avg{print}' ${perm_summary_stats_file} | wc -l)
    binding_diff_weaken_abs_avg_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)
    binding_diff_weaken_abs_avg_pval=$(echo -e "${binding_diff_weaken_abs_avg_num}\t${binding_diff_weaken_abs_avg_denom}" | awk -F"\t" '{print $1/$2}')
    
    perm_binding_diff_weaken_abs_avg_mean=$(awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_binding_diff_weaken_abs_avg_sd=$(awk '{sum += $6; sumsq += ($6)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})
    
    var=$(awk 'BEGIN{ print "'${perm_binding_diff_weaken_abs_avg_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
        binding_diff_weaken_abs_avg_z_score="NA"
    else
        binding_diff_weaken_abs_avg_z_score=$( echo -e "${real_binding_diff_weaken_abs_avg}\t${perm_binding_diff_weaken_abs_avg_mean})\t${perm_binding_diff_weaken_abs_avg_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi


    echo -e "${tf_name}\t${real_count}\t${count_pval}\t${perm_count_mean}\t${perm_count_sd}\t${count_z_score}\t${real_binding_diff_abs_avg}\t${binding_diff_abs_avg_pval}\t${perm_binding_diff_abs_avg_mean}\t${perm_binding_diff_abs_avg_sd}\t${binding_diff_abs_avg_z_score}\t${real_count_strengthen}\t${count_strengthen_pval}\t${perm_count_strengthen_mean}\t${perm_count_strengthen_sd}\t${count_strengthen_z_score}\t${real_binding_diff_strengthen_abs_avg}\t${binding_diff_strengthen_abs_avg_pval}\t${perm_binding_diff_strengthen_abs_avg_mean}\t${perm_binding_diff_strengthen_abs_avg_sd}\t${binding_diff_strengthen_abs_avg_z_score}\t${real_count_weaken}\t${count_weaken_pval}\t${perm_count_weaken_mean}\t${perm_count_weaken_sd}\t${count_weaken_z_score}\t${real_binding_diff_weaken_abs_avg}\t${binding_diff_weaken_abs_avg_pval}\t${perm_binding_diff_weaken_abs_avg_mean}\t${perm_binding_diff_weaken_abs_avg_sd}\t${binding_diff_weaken_abs_avg_z_score}" >> ${outFile}

done<${memeFileNamesFile}
