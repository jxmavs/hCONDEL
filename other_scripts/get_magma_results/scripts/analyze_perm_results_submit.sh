header=$1
analysis_ref_file=$2
script_dir=$3
temp_dir=$4

rm -f ${temp_dir}/${header}_analyze_perm_results.txt 
while read -a array
do
    analysis_run_name=${array[0]}
    column_metadata_file=${array[1]}
    empirical_stats_file=${array[2]}
    perm_stats_file=${array[3]}
    perm_analysis_out_file=${array[4]}

    echo -e "bash ${script_dir}/analyze_perm_results.sh '${analysis_run_name}' '${column_metadata_file}' '${empirical_stats_file}' '${perm_stats_file}' '${perm_analysis_out_file}'" >> ${temp_dir}/${header}_analyze_perm_results.txt
done<${analysis_ref_file}

#run everything

totalMem=2G
totalNumTasks=$(wc -l ${temp_dir}/${header}_analyze_perm_results.txt  | awk '{print $1}')   
totalTime=3:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${script_dir}/runTasksEval.sh > ${temp_dir}/${header}_analyze_perm_results.sh
rm -f ${temp_dir}/${header}_analyze_perm_results.log
qsub -o ${temp_dir}/${header}_analyze_perm_results.log ${temp_dir}/${header}_analyze_perm_results.sh ${temp_dir}/${header}_analyze_perm_results.txt

