file_names_file=$1
#all_subset_file_names_file contains all the file names of the subset files to combine
all_subset_file_names_file=$2
suffix=$3
#split_out_path is where the split files are sitting
split_out_path=$4
#out_file_path is where to store the combined file
out_file_path=$5
header_name=$6
temp_dir=$7

mkdir -p ${temp_dir}


rm -f ${temp_dir}/cat_merge_all_${header_name}.txt
#submit jobs
while read -a array
do

    file_name="${array[0]}"
    file_names_file_iter=${all_subset_file_names_file}.${file_name}.subset
    
    #get the subset of files corresponding to file_name    
    grep ${file_name} ${all_subset_file_names_file} | awk -v suffix=${suffix} -v split_out_path=${split_out_path} '{print split_out_path"/"$1"_"suffix}' > ${file_names_file_iter}
    
    output_file_name="${out_file_path}/${file_name}_${suffix}"
    echo "bash cat_merge.sh ${file_names_file_iter} ${output_file_name}" >> ${temp_dir}/cat_merge_all_${header_name}.txt

done<${file_names_file}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/cat_merge_all_${header_name}.txt  | awk '{print $1}')
totalTime=200:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/cat_merge_all_${header_name}.sh
rm -f ${temp_dir}/cat_merge_all_${header_name}.log
qsub -o ${temp_dir}/cat_merge_all_${header_name}.log  ${temp_dir}/cat_merge_all_${header_name}.sh ${temp_dir}/cat_merge_all_${header_name}.txt
