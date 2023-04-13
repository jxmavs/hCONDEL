#takes a list of files and gzips them

file_names_file=$1
header_name=$2
temp_dir=$3

mkdir -p ${temp_dir}

rm -f ${temp_dir}/gzip_submit_${header_name}.txt

script_dir=/cluster_path/ape_project/deletions_project/plasmid_seq/scripts

while read -a array
do

    file_path="${array[0]}"

    echo "bash ${script_dir}/gzip_iter.sh ${file_path}" >> ${temp_dir}/gzip_submit_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=20G
totalNumTasks=$(wc -l ${temp_dir}/gzip_submit_${header_name}.txt  | awk '{print $1}')
totalTime=30:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/gzip_submit_${header_name}.sh
rm -f ${temp_dir}/gzip_submit_${header_name}.log 
qsub -o ${temp_dir}/gzip_submit_${header_name}.log  ${temp_dir}/gzip_submit_${header_name}.sh ${temp_dir}/gzip_submit_${header_name}.txt
