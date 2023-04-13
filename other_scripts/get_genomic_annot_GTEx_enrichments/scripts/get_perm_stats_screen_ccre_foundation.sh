type_bed_file=$1
posFile=$2
out_path=$3
permFileHeader=$4

python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin
scriptdir=/cluster_path/ape_project/deletions_project/roadmap/scripts

del_dir=/cluster_path/ape_project/deletions_project

#intersect with screen_ccre
bedtools intersect -sorted -a ${posFile} -b ${type_bed_file} -loj > ${out_path}/${permFileHeader}_screen_ccre_intersected.txt

${python_dir}/python ${del_dir}/get_random_intersected_label.py ${out_path}/${permFileHeader}_screen_ccre_intersected.txt ${out_path}/${permFileHeader}.screen_ccre_info.txt

