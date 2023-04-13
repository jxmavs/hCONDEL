### get stats of permutations ###
tf_db_name=$1
fimo_out_path=$2
type_header_id=$3
tf_name=$4

#rm -f ${fimo_out_path}/perm_final/${type_header_id}_${tf_name}_summary_stats.txt
rm -f /cluster_path/ape_project/deletions_project/tf_enrichment/output/${type_header_id}_${tf_name}_summary_stats.txt
while read -a array
do
	perm_iter_file=${array[0]}
	#get count
	count=$(awk '($2!=0) {print}' ${perm_iter_file} | wc -l | awk '{print $1}')
	
	#get average binding score
	#NEED TO MAKE SURE NO final new line, otherwise it will print 0 REMEMBER TO CHECK OTHER FILES BEFORE?
	binding_diff_abs_sum=$(awk '(($2!=0) && ($15!="NA")) {print $15}' ${perm_iter_file} | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
	total_num_lines=$(wc -l ${perm_iter_file} | awk '{print $1}')			
	binding_diff_abs_avg=$(echo -e "${binding_diff_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')

    #get counts of binding site strengthening
    count_strengthen=$(awk '($2!=0) && ($15>0) {print}' ${perm_iter_file} | wc -l | awk '{print $1}')

	#get average binding score for strengthened sites
    binding_diff_strengthen_abs_sum=$(awk '(($2!=0) && ($15!="NA") && ($15>0)) {print $15}' ${perm_iter_file} | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
    binding_diff_strengthen_abs_avg=$(echo -e "${binding_diff_strengthen_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')
	
	#get counts of binding site weakening
	count_weaken=$(awk '($2!=0) && ($15<0) {print}' ${perm_iter_file} | wc -l | awk '{print $1}')	

    #get average binding score for weakened sites
    binding_diff_weaken_abs_sum=$(awk '(($2!=0) && ($15!="NA") && ($15<0)) {print $15}' ${perm_iter_file} | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1)}' | awk '{ sum += $1; n++ } END { if(n>0){print sum} else{print 0} }')
    binding_diff_weaken_abs_avg=$(echo -e "${binding_diff_weaken_abs_sum}\t${total_num_lines}" | awk -F"\t" '{print $1/$2}')
	
	#echo -e "${count}\t${binding_diff_abs_avg}\t${count_strengthen}\t${count_weaken}" >> ${fimo_out_path}/perm_final/${type_header_id}_${tf_name}_summary_stats.txt
	echo -e "${count}\t${binding_diff_abs_avg}\t${count_strengthen}\t${binding_diff_strengthen_abs_avg}\t${count_weaken}\t${binding_diff_weaken_abs_avg}" >> /cluster_path/ape_project/deletions_project/tf_enrichment/output/${type_header_id}_${tf_name}_summary_stats.txt

done<${fimo_out_path}/perm_final/${tf_name}_file_list.txt
