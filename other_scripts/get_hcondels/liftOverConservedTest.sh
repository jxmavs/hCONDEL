#!/bin/sh
#$ -cwd  
#$ -l mfree=7G
#$ -q long
#$ -V
#$ -pe smp 2
#$ -binding linear:2
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib



#te liftOver chain filminmatch=$1
minmatch=$1
phastconsdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/cons/run.phast

coord_dir=/cluster_path/ape_project/deletions_project/deletion_coordinates

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

del_dir=/cluster_path/ape_project/deletions_project

#awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$2"\t"$3}' ${phastconsdir}/mostConserved.bed > ${coord_dir}/mostConserved_meta.bed

#gunzip -c ${chimp_human_pairwise_path}/panTro4.hg38.all.chain.gz > ${chimp_human_pairwise_path}/panTro4.hg38.all.chain

#netChainSubset ${chimp_human_pairwise_path}/panTro4.hg38.net ${chimp_human_pairwise_path}/panTro4.hg38.all.chain ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain

#gzip ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain


#liftOver -minMatch=${minmatch} ${coord_dir}/mostConserved_meta.bed ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${coord_dir}/mostConservedMapped_new_meta_${minmatch}.txt ${coord_dir}/mostConservedNotMapped_new_meta_${minmatch}.txt

#liftOver ${phastconsdir}/mostConserved.bed ${del_dir}/chain_files/panTro4.hg38.all.chain.gz ${del_dir}/test/mostConservedMapped.txt ${del_dir}/test/mostConservedNotMapped.txt &

#liftOver ${phastconsdir}/mostConserved.bed ${chimp_human_pairwise_path}/panTro4.hg38.all.chain.gz ${del_dir}/test/mostConservedMapped_new.txt ${del_dir}/test/mostConservedNotMapped_new.txt &
parallel --link liftOver ::: -minMatch=${minmatch} ::: ${coord_dir}/mostConserved_meta.bed  ::: ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${del_dir}/chain_files/panTro4ToHg38.over.chain.gz ::: ${coord_dir}/mostConservedMapped_new_meta_${minmatch}.txt ${coord_dir}/mostConservedMapped_old_meta_${minmatch}.txt ::: ${coord_dir}/mostConservedNotMapped_new_meta_${minmatch}.txt ${coord_dir}/mostConservedNotMapped_old_meta_${minmatch}.txt

