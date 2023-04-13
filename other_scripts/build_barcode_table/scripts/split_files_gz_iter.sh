#https://www.gnu.org/software/coreutils/manual/html_node/split-invocation.html#split-invocation
file_name=$1
file_path=$2
split_out_path=$3
num_lines_split=$4

#the older versions of split do not have the filter option
split_bin_dir=/cluster_path/bin/coreutils_8.30/bin
gunzip -c ${file_path}/${file_name}.fastq.gz | ${split_bin_dir}/split -d -a 4 -l ${num_lines_split} --filter='gzip > $FILE.fastq.gz' - ${split_out_path}/${file_name}_
