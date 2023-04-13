tf_db_name=$1
fimo_out_path=$2
type_header_id=$3
tf_name=$4
num_lines=$5

rm -fr ${fimo_out_path}/perm_final/${tf_name}
mkdir -p ${fimo_out_path}/perm_final/${tf_name}
rm -f ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt

awk 'NR>1 {print}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.txt > ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt

#sort to reorder - added 11/7/21
awk '{print $1}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt | awk -F"_" '{print $NF}' > ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.perm_ids.txt
awk '{print $1}' ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt | awk -F"_" '{for (i = 1; i <= NF-2; i++) {printf $i"_"}; printf $(NF-1)"\n"}' > ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.orig_ids.txt

paste <(cat ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.perm_ids.txt) <(cat ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.orig_ids.txt) <(cat ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt) | sort -k1,1n -k2,2 | cut -f3- > ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.temp.txt

mv ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.temp.txt ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt
rm -f ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.perm_ids.txt
rm -f ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.orig_ids.txt 

split -d -l ${num_lines} -a 4 ${fimo_out_path}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.header_removed.txt ${fimo_out_path}/perm_final/${tf_name}/${type_header_id}_${tf_db_name}_fimo_cleaned_max_diff_output_final.${tf_name}.
ls -l ${fimo_out_path}/perm_final/${tf_name} | awk -v file_path=${fimo_out_path}/perm_final/${tf_name} 'NR>1{print file_path"/"$NF}' > ${fimo_out_path}/perm_final/${tf_name}_file_list.txt

