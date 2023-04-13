run_name=$1
ref_file=$2

temp_run_name=${run_name}_got

script_dir=/cluster_path/ape_project/deletions_project/10X/scripts

cd ${run_path}

temp_dir=/cluster_path_temp/10X
mkdir -p ${temp_dir}

rm -f ${temp_dir}/${temp_run_name}_all_runs.txt
#submit jobs
while read -a array
do
	config=${array[0]}
	fastqR1=${array[1]}
	fastqR2=${array[2]}
	whitelist=${array[3]}
	outdir=${array[4]}
	log=${array[5]}
	sample=${array[6]}
	
	echo "bash ${script_dir}/GoT_iter.sh ${config} ${fastqR1} ${fastqR2} ${whitelist} ${outdir} ${log} ${sample}" >> ${temp_dir}/${temp_run_name}_all_runs.txt

done<${ref_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
numCores=4
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/${temp_run_name}_all_runs.txt  | awk '{print $1}')
totalTime=100:00:00
sed "s/#$ -pe smp numCores/#$ -pe smp ${numCores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numCores}/g; s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/${temp_run_name}_submit_all.sh
rm -f ${temp_dir}/${temp_run_name}_submit_all.log
qsub -o ${temp_dir}/${temp_run_name}_submit_all.log ${temp_dir}/${temp_run_name}_submit_all.sh ${temp_dir}/${temp_run_name}_all_runs.txt
