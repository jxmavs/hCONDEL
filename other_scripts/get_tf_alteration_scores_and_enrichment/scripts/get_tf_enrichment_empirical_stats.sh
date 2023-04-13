### get stats of permutations ###
tf_db_name=$1
fimo_out_path=$2
memeFileNamesFile=$3
type_header_id=$4


### get stats of real data ###
#note to skip header line here
#rm -f ${fimo_out_path}/${type_header_id}_summary_stats.txt
rm -f /cluster_path/ape_project/deletions_project/tf_enrichment/output/${type_header_id}_summary_stats.txt

while read -a array
do
	motif_file=${array[0]}
	tf_name=$( echo ${motif_file} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' )
	echo ${tf_name}
	rm -f ${fimo_out_path}/${type_header_id}_${tf_name}_summary_stats.txt	
	#get count
	count=$(awk '(NR>1) && ($2!=0) {print}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | wc -l | awk '{print $1}')
	
	#get average binding score
	#NEED TO MAKE SURE NO final new line, otherwise it will print 0 REMEMBER TO CHECK OTHER FILES BEFORE?
	binding_diff_abs_sum=$(awk '((NR>1) && ($2!=0) && ($15!="NA")) {print $15}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
	total_num_lines=$(wc -l ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | awk '{print $1}')			
	binding_diff_abs_avg=$(echo -e "${binding_diff_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')
	
	#get counts of binding site strengthening
    count_strengthen=$(awk '(NR>1) && ($2!=0) && ($15>0) {print}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | wc -l | awk '{print $1}')

    #get average binding score for strengthened sites
    binding_diff_strengthen_abs_sum=$(awk '((NR>1) && ($2!=0) && ($15!="NA") && ($15>0)) {print $15}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
    binding_diff_strengthen_abs_avg=$(echo -e "${binding_diff_strengthen_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')

    #get counts of binding site weakening
    count_weaken=$(awk '(NR>1) && ($2!=0) && ($15<0) {print}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | wc -l | awk '{print $1}')

    #get average binding score for weakened sites
    binding_diff_weaken_abs_sum=$(awk '((NR>1) && ($2!=0) && ($15!="NA") && ($15<0)) {print $15}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
    binding_diff_weaken_abs_avg=$(echo -e "${binding_diff_weaken_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')

	#echo -e "${tf_name}\t${count}\t${binding_diff_abs_avg}\t${count_strengthen}\t${count_weaken}" >> ${fimo_out_path}/${type_header_id}_summary_stats.txt

	echo -e "${tf_name}\t${count}\t${binding_diff_abs_avg}\t${count_strengthen}\t${binding_diff_strengthen_abs_avg}\t${count_weaken}\t${binding_diff_weaken_abs_avg}" >> /cluster_path/ape_project/deletions_project/tf_enrichment/output/${type_header_id}_summary_stats.txt
	
done<${memeFileNamesFile}


