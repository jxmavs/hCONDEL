#!/bin/sh
#$ -cwd 
#$ -l h_vmem=2G
#$ -q long
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7

LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.2.1/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib


#download roadmap data

inputdir=$1
outdir=$2
modType=$3
roadFileType=$4
IDListFile=$5

wrkdir=/cluster_path_temp/roadmap

rm -f ${wrkdir}/all_${roadFileType}_bigwig_decomp_${modType}.txt

while read -a array
do
	roadmapID=${array[0]}
    wiggleFileName=${inputdir}/${roadmapID}-${modType}.${roadFileType}.signal.bigwig
    outFileName=${outdir}/${roadmapID}-${modType}.${roadFileType}.signal.bedgraph   
    echo "bash decompress_bigwig_iter.sh ${wiggleFileName} ${outFileName}" >> ${wrkdir}/all_${roadFileType}_bigwig_decomp_${modType}.txt
	
done<${IDListFile}

#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=2G
totalNumTasks=$(wc -l ${wrkdir}/all_${roadFileType}_bigwig_decomp_${modType}.txt  | awk '{print $1}')   
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${wrkdir}/decompress_${roadFileType}_bigwig_${modType}_all.sh
rm -f ${wrkdir}/decompress_${roadFileType}_bigwig_${modType}_all.log
qsub -o ${wrkdir}/decompress_${roadFileType}_bigwig_${modType}_all.log ${wrkdir}/decompress_${roadFileType}_bigwig_${modType}_all.sh ${wrkdir}/all_${roadFileType}_bigwig_decomp_${modType}.txt
