#!/bin/sh
#$ -cwd 
#$ -l m_mem_free=5G
#$ -q long
#$ -V
#$ -j y
source /broad/software/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib

TNAME=$1
QNAME=$2
WRKDIR=$3
OUTDIR=$4
#this script merges all chains produced

#remove psl files from lastz
#rm -r ${WRKDIR}/psl

#combine all chain files into larger file, make sure to add chainPreNet!!!!
find ${WRKDIR}/chain -name "*.chain" | chainMergeSort -inputList=stdin | gzip -c > ${OUTDIR}/${TNAME}.${QNAME}.all.chain.gz

#remove original chain files
#rm -r ${WRKDIR}/chain

