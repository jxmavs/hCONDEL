file_names_file=$1
file_path=$2
split_out_path=$3
num_lines_split=$4
header_name=$5
temp_dir=$6
#############################################################################################################

mkdir -p ${temp_dir}

#how much do we trim off beginning/end, remove barcode and adapter sequences
#cutAMTEnd removes CGTCAAGCGGCCAGTT (16 bp) adapter sequence from end
#cutAMTBegin removes barcode (20bp barcode sequence) + adapter TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG (38 bp long) 

rm -f ${temp_dir}/split_files_all_${header_name}.txt
#submit jobs
while read -a array
do

    file_name="${array[0]}"
    echo "bash split_files_gz_iter.sh ${file_name} ${file_path} ${split_out_path} ${num_lines_split}" >> ${temp_dir}/split_files_all_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/split_files_all_${header_name}.txt  | awk '{print $1}')
totalTime=200:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/split_files_all_${header_name}.sh
rm -f ${temp_dir}/split_files_all_${header_name}.log
qsub -o ${temp_dir}/split_files_all_${header_name}.log  ${temp_dir}/split_files_all_${header_name}.sh ${temp_dir}/split_files_all_${header_name}.txt
