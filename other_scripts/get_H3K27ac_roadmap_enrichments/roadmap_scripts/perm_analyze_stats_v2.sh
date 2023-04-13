#combine all permutation files into individual cell file

modType=$1
cutOffVal=$2
real_stat_file=$3
roadmap_ID_file=$4
out_header=$5
numPerm=$6
outDir=$7
outFile=${outDir}/all_${out_header}_${modType}_cutoff_${cutOffVal}_final_perm_stats.txt

rm -f ${outFile}
while IFS=$'\t' read -a array
do
    echo ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt
    roadmapID=${array[0]}
    cellType=${array[1]}   
    
    
    ######### for p-values ##############
    #get real value
    real_avg_score_pval=$(awk -F"\t" -v roadmapID=${roadmapID} '$(NF-1)==roadmapID {print $1}' ${real_stat_file})
    avg_score_num_pval=$(awk -F"\t" -v real_avg_score_pval=${real_avg_score_pval} '$1>real_avg_score_pval{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    avg_score_denom_pval=$(awk '{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    avg_score_val_pval=$(echo -e "${avg_score_num_pval}\t${avg_score_denom_pval}" | awk -F"\t" '{print $1/$2}')    
    
    perm_avg_score_pval=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	perm_sd_score_pval=$(awk '{sum += $1; sumsq += ($1)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	 
    #get real value
    real_prop_int_pval=$(awk -F"\t" -v roadmapID=${roadmapID} '$(NF-1)==roadmapID {print $2}' ${real_stat_file})
    prop_int_num_pval=$(awk -F"\t" -v real_prop_int_pval=${real_prop_int_pval} '$2>real_prop_int_pval{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l )
    prop_int_denom_pval=$(awk '{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    prop_int_val_pval=$(echo -e "${prop_int_num_pval}\t${prop_int_denom_pval}" | awk -F"\t" '{print $1/$2}')    

    perm_avg_prop_int_pval=$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	perm_sd_prop_int_pval=$(awk '{sum += $2; sumsq += ($2)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
    
	########## for fc #######################

    real_avg_score_fc=$(awk -F"\t" -v roadmapID=${roadmapID} '$(NF-1)==roadmapID {print $4}' ${real_stat_file})
    avg_score_num_fc=$(awk -F"\t" -v real_avg_score_fc=${real_avg_score_fc} '$4>real_avg_score_fc{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    avg_score_denom_fc=$(awk '{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    avg_score_val_fc=$(echo -e "${avg_score_num_fc}\t${avg_score_denom_fc}" | awk -F"\t" '{print $1/$2}')    
    
    perm_avg_score_fc=$(awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	perm_sd_score_fc=$(awk '{sum += $4; sumsq += ($4)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	
    #get real value
    real_prop_int_fc=$(awk -F"\t" -v roadmapID=${roadmapID} '$(NF-1)==roadmapID {print $5}' ${real_stat_file})
    prop_int_num_fc=$(awk -F"\t" -v real_prop_int_fc=${real_prop_int_fc} '$5>real_prop_int_fc{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l )
    prop_int_denom_fc=$(awk '{print}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt | wc -l)
    prop_int_val_fc=$(echo -e "${prop_int_num_fc}\t${prop_int_denom_fc}" | awk -F"\t" '{print $1/$2}')    

    perm_avg_prop_int_fc=$(awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)
	perm_sd_prop_int_fc=$(awk '{sum += $5; sumsq += ($5)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt)

    #echo -e "${roadmapID}\t${cellType}\t${avg_score_val_pval}\t${prop_int_val_pval}\t${avg_score_val_fc}\t${prop_int_val_fc}\t${real_avg_score_pval}\t${real_prop_int_pval}\t${real_avg_score_fc}\t${real_prop_int_fc}\t${perm_avg_score_pval}\t${perm_avg_prop_int_pval}\t${perm_avg_score_fc}\t${perm_avg_prop_int_fc}" 
    echo -e "${roadmapID}\t${cellType}\t${avg_score_val_pval}\t${prop_int_val_pval}\t${avg_score_val_fc}\t${prop_int_val_fc}\t${real_avg_score_pval}\t${real_prop_int_pval}\t${real_avg_score_fc}\t${real_prop_int_fc}\t${perm_avg_score_pval}\t${perm_avg_prop_int_pval}\t${perm_avg_score_fc}\t${perm_avg_prop_int_fc}\t${perm_sd_score_pval}\t${perm_sd_prop_int_pval}\t${perm_sd_score_fc}\t${perm_sd_prop_int_fc}"  >> ${outFile}
done<${roadmap_ID_file}
