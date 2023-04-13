type_iter=$1
type_perm_stats_file_list_file=$2
perm_type_results_dir=$3

type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_random
rm -f ${perm_type_results_dir}/${type_iter}/${type_header_id}_score_summary_stats.txt

while read -a array
do
	perm_iter_file=${array[0]}
	if [[ ${type_iter} =~ "gtex_v8" ]]; then
		perm_score_avg=$(awk '$NF>-1{ sum += $NF; n++ } END { if (n > 0) print sum / n; }' ${perm_iter_file})
	else
		perm_score_avg=$(awk '{ sum += $NF; n++ } END { if (n > 0) print sum / n; }' ${perm_iter_file})
	fi
	
	echo -e "${perm_score_avg}" >> ${perm_type_results_dir}/${type_iter}/${type_header_id}_score_summary_stats.txt

done<${type_perm_stats_file_list_file}
