type_bed_file=$1
posFile=$2
out_path=$3
permFileHeader=$4

script_dir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

dist_cutoff_interest_file="/cluster_path/ape_project/deletions_project/type/type_dist_cutoff_interest_file.txt"

#get closest feature, parse to see if it passes distance cutoffs in ${dist_cutoff_interest_file} 
bedtools closest -d -a ${posFile} -b ${type_bed_file} | awk '{print $4"\t"$NF}' > "${out_path}/${permFileHeader}.sample_info.temp.txt"
python ${script_dir}/format_bedtools_closest_dist_output.py "${out_path}/${permFileHeader}.sample_info.temp.txt" "${dist_cutoff_interest_file}" "${out_path}/${permFileHeader}.sample_info.txt"
rm -f "${out_path}/${permFileHeader}.sample_info.temp.txt"
