#!/bin/bash
#$ -cwd 
#$ -l h_vmem=5G
#$ -l h_rt=30:00:00
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
reuse -q Perl-5.20
reuse -q BLAST+
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/cluster_path/bin/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/mysql_5.6.20/lib

type_ref_file=$1
perm_type_results_dir=$2
script_dir=$3

while read -a array
do
    type_iter=${array[0]}
    type_bed_file=${array[1]}
    inputPermFileList=${array[4]}

    out_path=${perm_type_results_dir}/${type_iter}
    rm -fr ${out_path}
    mkdir -p ${out_path}

    out_header=${type_iter}_enrichment_hcondel_denylist_filtered_random

    temp_dir=/cluster_path_temp/${type_iter}_perm_stats
    rm -fr ${temp_dir}
    mkdir -p ${temp_dir}

    #submit
    bash ${script_dir}/get_perm_stats_submit_v2.sh ${type_bed_file} ${out_path} ${out_header} ${inputPermFileList} ${temp_dir} ${type_iter}

done<${type_ref_file}

