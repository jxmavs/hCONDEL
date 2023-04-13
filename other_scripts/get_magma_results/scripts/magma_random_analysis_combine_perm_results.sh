random_analysis_header=$1
gwas_analysis_header=$2
gwas_analyses_of_interest_metadata_file=$3
batch_num_file=$4
magma_random_analyses_type_dir=$5

rm -f ${magma_random_analyses_type_dir}/${gwas_analysis_header}/${random_analysis_header}_${gwas_analysis_header}_all_perm_stats.txt

while IFS=$'\t' read -a array
do
	
	batch_num="${array[0]}"
	analysis_header=${magma_random_analyses_type_dir}/${gwas_analysis_header}/${random_analysis_header}_${gwas_analysis_header}_batch_${batch_num}

	grep -v "#" ${analysis_header}_gs_condition_cons_prop.gsa.out | awk '$1=="Overlap" && $3=="2"{print $5"\t"$6"\t"$7"\t"$8}' >> ${magma_random_analyses_type_dir}/${gwas_analysis_header}/${random_analysis_header}_${gwas_analysis_header}_all_perm_stats.txt

done<${batch_num_file}

