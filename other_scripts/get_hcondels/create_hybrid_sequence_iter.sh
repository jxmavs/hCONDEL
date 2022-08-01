#!/bin/sh
#$ -cwd 
#$ -l h_vmem=30G
#$ -q long
#$ -V
#$ -j y

source /broad/software/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.2.1/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib


chimp_insert_coord=$1
chimp_insert_fasta=$2
human_ref_fasta_file=$3
human_chr_sizes_file=$4
chimp_human_hybrid_out_fasta=$5
chimp_human_hybrid_out_meta=$6
#combine conserved with deleted sequence
#old code without taking into account reverse complement
#python createHybridSequence.py ${chimp_insert_coord} ${chimp_insert_fasta} ${human_ref_fasta_file} ${human_chr_sizes_file} ${chimp_human_hybrid_out_fasta} ${chimp_human_hybrid_out_meta}


python createHybridSequenceWithRevComp.py ${chimp_insert_coord} ${chimp_insert_fasta} ${human_ref_fasta_file} ${human_chr_sizes_file} ${chimp_human_hybrid_out_fasta} ${chimp_human_hybrid_out_meta}

