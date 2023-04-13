scriptdir=$1
out_path=$2
pairwise_comp_file=$3
tf_enrichment_files_dir=$4
type_file=$5
tf_db_name=$6
fimo_out_path=$7

del_dir=/cluster_path/ape_project/deletions_project

rm -f ${del_dir}/tf_enrichment_remove_files_all_runs.txt

while read -a array1
do
	type_header_id=${array1[0]}
	seqFileNamesFile=${array1[1]}
	memeFileNamesFile=${array1[2]}

	while read -a array2
	do
		motif_file=${array2[0]}
		tf_name=$( echo ${motif_file} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' )
			
		#out_header=${cell_header_id}_del_encode_annotation_${tf_db_name}
		#file_out_header=encode_annote_${tf_db_name}
		#bash ${scriptdir}/fimo_combine_all.sh ${out_header} ${file_out_header} ${del_dir} ${pairwise_comp_file}
		
		out_header=${type_header_id}_${tf_db_name}
		file_out_header=${tf_db_name}
		echo -e "bash ${scriptdir}/tf_enrichment_remove_files_iter.sh ${out_header} ${file_out_header} ${fimo_out_path} ${pairwise_comp_file} '${tf_name}' ${seqFileNamesFile}" >> ${del_dir}/tf_enrichment_remove_files_all_runs.txt
		
		#sort by date to double check correct files
		#ls -t --full-time  ${del_dir}/${out_header}_* | tail
		
	done<${memeFileNamesFile}

done<${type_file}


del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${del_dir}/tf_enrichment_remove_files_all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${del_dir}/tf_enrichment_remove_files_all_runs.sh
rm -f ${del_dir}/tf_enrichment_remove_files_all_runs.log 
qsub -o ${del_dir}/tf_enrichment_remove_files_all_runs.log ${del_dir}/tf_enrichment_remove_files_all_runs.sh ${del_dir}/tf_enrichment_remove_files_all_runs.txt



