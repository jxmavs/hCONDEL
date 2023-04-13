type=$1
ls -l /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/*pval.signal.bigwig | awk '{print $NF}' > /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_initial.txt
rm -f /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_temp.txt
while read -a array
do
    file_name=${array[0]}
    if [ -s ${file_name} ]
    then
        echo ${file_name} >> /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_temp.txt
    else
        echo ${file_name}
    fi
done</cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_initial.txt

temp_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_temp.txt
paste <(awk -F"/" '{print $NF}'  ${temp_file} | awk -F"-" '{print $1}' ) <( awk '{print}' ${temp_file} ) > /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
head /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt

