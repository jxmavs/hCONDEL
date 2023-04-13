#!/bin/bash
#$ -cwd 
#$ -l h_vmem=totalMem
#$ -l h_rt=totalTime
#$ -t 1-totalNumTasks
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
reuse -q Perl-5.20
reuse -q BLAST+
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/cluster_path/bin/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/mysql_5.6.20/lib

#this script submits a task array and evaluates each of the tasks from a text file
tasklist=$1
runNum=${SGE_TASK_ID}

#run the command from the file
eval $(awk -v runNum=${runNum} 'NR==runNum{print}' ${tasklist})
