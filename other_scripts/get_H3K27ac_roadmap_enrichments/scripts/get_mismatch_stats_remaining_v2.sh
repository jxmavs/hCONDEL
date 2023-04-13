header_name=$1
out_path=$2
seq_table_file_path=$3
remaining_file_path=$4

rm -f ${remaining_file_path}

while read -a array 
do
    seq_id="${array[3]}"

    if [[ ! -f ${out_path}/${header_name}.human_coord.${seq_id}.blast.stats.parsed ]];
    then
        echo ${seq_id} >> ${remaining_file_path}
    fi

done<${seq_table_file_path}


