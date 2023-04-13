read_length=$1
fragment_length=$2
read_one_file=$3
read_two_file=$4
output_dir=$5
output_prefix=$6

flash --compress-prog=gzip -r ${read_length} -f ${fragment_length} -d ${output_dir} -o ${output_prefix} ${read_one_file} ${read_two_file}
