#submit jobs
header_name=$1
file_names_file=$2
read_length=$3
fragment_length=$4
output_dir=$5
temp_dir=$6

mkdir -p ${temp_dir}

rm -f ${temp_dir}/flash_merge_submit_${header_name}.txt

script_dir=/cluster_path/ape_project/deletions_project/plasmid_seq/scripts

while read -a array
do

    file_header="${array[0]}"
    read_one_file="${array[1]}"
    read_two_file="${array[2]}"

    echo "bash ${script_dir}/flash_merge_gz.sh ${read_length} ${fragment_length} ${read_one_file} ${read_two_file} ${output_dir} ${file_header}" >> ${temp_dir}/flash_merge_submit_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=50G
totalNumTasks=$(wc -l ${temp_dir}/flash_merge_submit_${header_name}.txt  | awk '{print $1}')
totalTime=200:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/flash_merge_submit_${header_name}.sh
rm -f ${temp_dir}/flash_merge_submit_${header_name}.log 
qsub -o ${temp_dir}/flash_merge_submit_${header_name}.log ${temp_dir}/flash_merge_submit_${header_name}.sh ${temp_dir}/flash_merge_submit_${header_name}.txt
