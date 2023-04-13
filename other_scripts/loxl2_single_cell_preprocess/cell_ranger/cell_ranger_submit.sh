run_name=$1
run_path=$2
ref_file=$3

temp_run_name=${run_name}_cellranger

script_dir=/cluster_path/ape_project/deletions_project/10X/scripts

cd ${run_path}

temp_dir=/cluster_path_temp/10X
mkdir -p ${temp_dir}

rm -f ${temp_dir}/${temp_run_name}_all_runs.txt
#submit jobs
while read -a array
do
	id=${array[0]}
	transcriptome=${array[1]}
	fastqs=${array[2]}
	sample=${array[3]}
	expect_cells=${array[4]}
	localcores=${array[5]}
	localmem=${array[6]}
	
	echo "bash ${script_dir}/cell_ranger_iter.sh ${id} ${transcriptome} ${fastqs} ${sample} ${expect_cells} ${localcores} ${localmem}" >> ${temp_dir}/${temp_run_name}_all_runs.txt

done<${ref_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
numCores=2
totalMem=20G
totalNumTasks=$(wc -l ${temp_dir}/${temp_run_name}_all_runs.txt  | awk '{print $1}')
totalTime=30:00:00
sed "s/#$ -pe smp numCores/#$ -pe smp ${numCores}/g; s/#$ -binding linear:numCores/#$ -binding linear:${numCores}/g; s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksParallel.sh > ${temp_dir}/${temp_run_name}_submit_all.sh
rm -f ${temp_dir}/${temp_run_name}_submit_all.log
qsub -o ${temp_dir}/${temp_run_name}_submit_all.log ${temp_dir}/${temp_run_name}_submit_all.sh ${temp_dir}/${temp_run_name}_all_runs.txt
