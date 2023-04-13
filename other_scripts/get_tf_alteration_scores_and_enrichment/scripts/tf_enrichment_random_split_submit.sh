scriptdir=$1
out_path=$2
tf_db_name=$3
fimo_out_path=$4
memeFileNamesFile=$5
type_header_id=$6

rm -f ${out_path}/tf_enrichment_random_split_all_runs.txt

while read -a array
do
	motif_file=${array[0]}
	tf_name=$( echo ${motif_file} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' )
	num_lines=$(wc -l /cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt | awk '{print $1}')
		
	echo -e "bash ${scriptdir}/tf_enrichment_random_split_iter.sh ${tf_db_name} ${fimo_out_path} ${type_header_id} '${tf_name}' ${num_lines}" >> ${out_path}/tf_enrichment_random_split_all_runs.txt
	
done<${memeFileNamesFile}



del_dir=/cluster_path/ape_project/deletions_project
totalMem=20G
totalNumTasks=$(wc -l ${out_path}/tf_enrichment_random_split_all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${out_path}/tf_enrichment_random_split_all_runs.sh
rm -f ${out_path}/tf_enrichment_random_split_all_runs.log 
qsub -o ${out_path}/tf_enrichment_random_split_all_runs.log ${out_path}/tf_enrichment_random_split_all_runs.sh ${out_path}/tf_enrichment_random_split_all_runs.txt


