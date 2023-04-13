#v5_altered splits into running each FIMO TF seperately 

out_header=$1
pairwise_comp_file=$2
lookAtCenteredOnly=$3
cutOffVal=$4
memeFileName=$5
out_path=$6

seqFileNamesFile=$7
memeFileNamesFile=$8

scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

rm -f ${out_path}/${out_header}_fimo_all_runs.txt
while read -a array
do
	input_file=${array[0]}
	while read -a array2
	do
		motif_file=${array2[0]}
		echo "bash ${scriptdir}/fimo_run_iter_v6_altered.sh ${input_file} ${out_header} '${motif_file}' ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${memeFileName} ${out_path}" >> ${out_path}/${out_header}_fimo_all_runs.txt
	done<${memeFileNamesFile}
done<${seqFileNamesFile}

#now output


del_dir=/cluster_path/ape_project/deletions_project
totalMem=20G
totalNumTasks=$(wc -l ${out_path}/${out_header}_fimo_all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${out_path}/${out_header}_fimo_all_runs.sh
rm -f ${out_path}/${out_header}_fimo_all_runs.log 
qsub -o ${out_path}/${out_header}_fimo_all_runs.log ${out_path}/${out_header}_fimo_all_runs.sh ${out_path}/${out_header}_fimo_all_runs.txt



