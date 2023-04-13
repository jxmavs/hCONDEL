#!/bin/sh
#$ -cwd 
#$ -l h_vmem=5G
#$ -q short
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib


file_name=$1
file_path=$2
temp_path=$3
out_path=$4
bowtie_ind=$5
file_suffix=$6
numCores=$7

newFileName=${file_name}_${file_suffix}

python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin

#trim the fastq files

#cutadapt -u -${cutAMT} -o ${temp_path}/${newFileName}_trimmed.fastq ${file_path}/${file_name}.fastq.gz
#cut off first bp


#gunzip -c ${file_path}/${file_name}.fastq.gz > ${temp_path}/${file_name}.fastq
#add on barcode info also to fastq
#gunzip -c ${file_path}/${file_name}.fastq.gz | awk -v cutAMTBegin=${cutAMTBegin} -v cutAMTEnd=${cutAMTEnd} '{if(NR%2!=0){print $1}else{print substr($1,cutAMTBegin+1, length($1)-cutAMTBegin-cutAMTEnd)}}' | gzip > ${temp_path}/${newFileName}_trimmed.fastq.gz

#add on a barcode, and trim, use leviathan distance to id barcode
#since the deletion product is 274 bp long, I can't have any deleted sequence that removes 124bp otherwise the sequence will be incorrect, since I take 20 bp from end of sequence as barcode

${python_dir}/python rmAdaptersV2MPRADel.py ${file_path}/${file_name}.fastq.gz -o ${out_path}/${newFileName}_trimmed.fastq.gz -d ${out_path}/${newFileName}_trimmed.discarded.fastq.gz -barxN -Ladp_trim_bp

#bowtie -v ${mismatchAMT} -m ${multimatchAMT} -S -q ${bowtie_ind} ${temp_path}/${newFileName}_trimmed.fastq | samtools view -bS > ${temp_path}/${newFileName}_trimmed.bam

#bowtie -v ${mismatchAMT} -m ${multimatchAMT} -S -q ${bowtie_ind} -1 ${temp_path}/${newFileName}_trimmed_1.fastq -2 ${temp_path}/${newFileName}_trimmed_2.fastq > ${temp_path}/${newFileName}_trimmed.sam

#bowtie -m ${multimatchAMT} -S -q ${bowtie_ind} -1 ${temp_path}/${newFileName}_trimmed_1.fastq -2 ${temp_path}/${newFileName}_trimmed_2.fastq > ${temp_path}/${newFileName}_trimmed.sam

#bowtie -m ${multimatchAMT} -S -q ${bowtie_ind} ${temp_path}/${newFileName}_trimmed.fastq >  ${temp_path}/${newFileName}_trimmed_bowtie.sam

#bowtie2 -x ${bowtie_ind} --very-sensitive -p ${numCores} -U ${temp_path}/${newFileName}_trimmed.fastq.gz --al-gz ${temp_path}/${newFileName}_trimmed_bowtie2.aligned.sam --un-gz ${temp_path}/${newFileName}_trimmed_bowtie2.unaligned.sam 2>${temp_path}/${newFileName}_stderr.txt > /dev/null

bowtie2 -x ${bowtie_ind} --very-sensitive -p ${numCores} -U ${out_path}/${newFileName}_trimmed.fastq.gz --al-gz ${out_path}/${newFileName}_trimmed_bowtie2.aligned.fastq.gz 2>${out_path}/${newFileName}_stderr.txt | gzip > ${out_path}/${newFileName}_trimmed_bowtie2.sam.gz

#seperate out aligned and unalinged reads 

gunzip -c ${out_path}/${newFileName}_trimmed_bowtie2.sam.gz | samtools view -b -F 4 - > ${out_path}/${newFileName}_trimmed_bowtie2.aligned.bam

gunzip -c ${out_path}/${newFileName}_trimmed_bowtie2.sam.gz | samtools view -b -f 4 - > ${out_path}/${newFileName}_trimmed_bowtie2.unaligned.bam

rm -f ${out_path}/${newFileName}_trimmed_bowtie2.sam.gz

#rm ${temp_path}/${newFileName}_*
#rm ${temp_path}/${file_name}_*

