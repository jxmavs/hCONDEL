#!/bin/sh
#$ -cwd 
#$ -pe smp 3
#$ -binding linear:3
#$ -l h_vmem=10G
#$ -q long
#$ -V
#$ -j y

source /broad/software/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.2.1/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib


sgdp_output_dir=$1
ref_file_path=$2
pop_sample_file_path=$3
pop_sample_id=$4
genome_file_id=$5
numcores=$6

fermi_program_path=/cluster_path/bin/fermikit/fermi.kit

${fermi_program_path}/bwa mem -t ${numcores} -x intractg ${ref_file_path} ${pop_sample_file_path} | gzip -1 > ${sgdp_output_dir}/${pop_sample_id}_${genome_file_id}.sam.gz
${fermi_program_path}/htsbox samsort -t ${numcores} -S ${sgdp_output_dir}/${pop_sample_id}_${genome_file_id}.sam.gz > ${sgdp_output_dir}/${pop_sample_id}_${genome_file_id}.srt.bam
rm ${sgdp_output_dir}/${pop_sample_id}_${genome_file_id}.sam.gz  
