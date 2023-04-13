type_bed_file=$1
out_path=$2
out_header=$3
inputPermFileList=$4
temp_dir=$5
type_iter=$6

script_dir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

permIter=1
rm -f ${temp_dir}/${out_header}_all_runs.txt 
while read -a array
do
	posFile="${array[0]}"
	permFileHeader=${out_header}_${permIter}
	permIter=$((${permIter} + 1)) 
	echo "bash '${script_dir}/get_perm_stats_${type_iter}.sh' '${type_bed_file}' '${posFile}' '${out_path}' '${permFileHeader}'" >> ${temp_dir}/${out_header}_all_runs.txt 
done<${inputPermFileList}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/${out_header}_all_runs.txt  | awk '{print $1}')   
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/${out_header}_getTypeScoresAll.sh
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/${out_header}_getTypeScoresAll.sh
rm -f ${temp_dir}/${out_header}_getTypeScoresAll.log
qsub -o ${temp_dir}/${out_header}_getTypeScoresAll.log ${temp_dir}/${out_header}_getTypeScoresAll.sh ${temp_dir}/${out_header}_all_runs.txt

