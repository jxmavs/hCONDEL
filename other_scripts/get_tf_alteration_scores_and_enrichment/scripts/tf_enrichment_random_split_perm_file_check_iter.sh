fimo_out_path=$1  
tf_name=$2

#rm -f ${fimo_out_path}/perm_final/${tf_name}_file_check.txt
rm -f /cluster_path/ape_project/deletions_project/tf_enrichment/output/${tf_name}_file_check.txt
while read -a array2
do
	perm_split_iter_file=${array2[0]}
	#wc -l ${perm_split_iter_file} >> ${fimo_out_path}/perm_final/${tf_name}_file_check.txt
	wc -l ${perm_split_iter_file} >> /cluster_path/ape_project/deletions_project/tf_enrichment/output/${tf_name}_file_check.txt
done<${fimo_out_path}/perm_final/${tf_name}_file_list.txt


