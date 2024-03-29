#!/bin/sh
#$ -cwd 
#$ -l h_vmem=5G
#$ -l h_rt=200:00:00
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib

file_names_file=$1
all_subset_file_names_file=$2
out_path=$3
header_name=$4
suffix=$5

rm -f ${out_path}/${header_name}_ehb_stats.txt

while read -a array
do
    file_header="${array[0]}"
    num_flash_merged=$( gunzip -c ${out_path}/${file_header}.fastq.gz | wc -l | awk '{print $1/4}' )
    num_flash_failed=$( gunzip -c ${out_path}/${file_header}.notCombined_1.fastq.gzip | wc -l | awk '{print $1/4}' )

    #file_name="${array[0]}"
    file_names_file_iter=${all_subset_file_names_file}.${file_header}.subset.header

    #get the subset of files corresponding to file_name    
    grep ${file_header} ${all_subset_file_names_file} | awk '{print $1}' > ${file_names_file_iter}
    #if no split files, then the file wasn't split 
    if [[ ! -s ${file_names_file_iter} ]]; 
    then
        #if the file is split, then sum up all the info from the split files
          
        newFileName=${file_header}_${suffix}

        #number of reads discraded due to not being able to locate the barcode
        num_barcode_adp_discarded=$( gunzip -c ${out_path}/${newFileName}.discarded.fastq.gz | wc -l | awk '{print $1/4}' )
     
        num_aligned=$( samtools view -c ${out_path}/${newFileName}_bowtie2.aligned.bam )
        num_unaligned=$( samtools view -c ${out_path}/${newFileName}_bowtie2.unaligned.bam )
    else
       
        num_aligned=0
        num_unaligned=0
        num_barcode_adp_discarded=0
         
        while read -a array
        do
            file_header_iter="${array[0]}"
            #if the file is split, then sum up all the info from the split files
         
            newFileName=${file_header_iter}_${suffix}

            #number of reads discraded due to not being able to locate the barcode
            num_barcode_adp_discarded_iter=$( gunzip -c ${out_path}/${newFileName}.discarded.fastq.gz | wc -l | awk '{print $1/4}' )
            
            num_aligned_iter=$( samtools view -c ${out_path}/${newFileName}_bowtie2.aligned.bam )
            num_unaligned_iter=$( samtools view -c ${out_path}/${newFileName}_bowtie2.unaligned.bam )
            
            num_aligned=$((num_aligned_iter+num_aligned))
            num_unaligned=$((num_unaligned_iter+num_unaligned))
            num_barcode_adp_discarded=$((num_barcode_adp_discarded_iter+num_barcode_adp_discarded))
            
        done<${file_names_file_iter}
    
    fi
    
    echo -e "${file_header}\t${num_flash_merged}\t${num_flash_failed}\t${num_barcode_adp_discarded}\t${num_aligned}\t${num_unaligned}" >> ${out_path}/${header_name}_ehb_stats.txt

done<${file_names_file}

#num_flash_merged+num_flash_failed is the initial number of reads
