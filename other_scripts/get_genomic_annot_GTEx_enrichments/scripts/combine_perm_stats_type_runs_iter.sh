type_iter=$1
perm_type_results_dir=$2

out_path=${perm_type_results_dir}/${type_iter}
ls -l ${out_path}/*${type_iter}_info.txt | awk '{print $NF}' | sort -k1,1 > ${perm_type_results_dir}/${type_iter}_perm_initial_results_file_list.txt

