#v5 is modified to change fimo_run_iter_v5.sh to fimo_run_iter_v6.sh, to allow for encode ids to be correctly matched
#obtain file list names
#ls -l /cluster_path_temp/del_encode/coord_split/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.* > ${out_path}/all_fimo_file_names.txt

out_header=$1
motif_file=$2
pairwise_comp_file=$3
lookAtCenteredOnly=$4
cutOffVal=$5
useEncodeAnnotation=$6
encodeFileName=$7
encodeIDMapFileName=$8
lookAtExpressedTFs=$9
expressedTFsFileName=${10}
memeFileName=${11}
out_path=${12}

scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

rm -f ${out_path}/${out_header}_fimo_all_runs.txt
while read -a array
do
    
    input_file=${array[0]}
    echo "bash ${scriptdir}/fimo_run_iter_v6.sh ${input_file} ${out_header} ${motif_file} ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${useEncodeAnnotation} ${encodeFileName} ${encodeIDMapFileName} ${lookAtExpressedTFs} ${expressedTFsFileName} ${memeFileName} ${out_path}" >> ${out_path}/${out_header}_fimo_all_runs.txt
done<${out_path}/all_fimo_file_names.txt

#now output


del_dir=/cluster_path/ape_project/deletions_project
totalMem=10G
totalNumTasks=$(wc -l ${out_path}/${out_header}_fimo_all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${out_path}/${out_header}_fimo_all_runs.sh
rm -f ${out_path}/${out_header}_fimo_all_runs.log 
qsub -o ${out_path}/${out_header}_fimo_all_runs.log ${out_path}/${out_header}_fimo_all_runs.sh ${out_path}/${out_header}_fimo_all_runs.txt



