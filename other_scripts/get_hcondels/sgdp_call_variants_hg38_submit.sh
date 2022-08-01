#!/bin/sh
#$ -cwd 
#$ -l h_vmem=20G
#$ -q long
#$ -V
#$ -j y

source /broad/software/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.2.1/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib

sgdp_bam_dir=$1
sgdp_vcf_dir=$2
ref_fasta_name=$3
sgdp_indiv_names_file=$4
genome_file_id=$5

fermi_program_path=/cluster_path/bin/fermikit/fermi.kit
temp_dir=/cluster_path_temp/sgdp/temp

#this script calls all short and large variants
#calls all short variants together with all bam files
qsub call_sgdp_small_variants.sh ${ref_fasta_name} ${sgdp_bam_dir} ${sgdp_vcf_dir} ${genome_file_id}

#call all short variants
#${fermi_program_path}/htsbox pileup -cuf ${ref_file_path} ${sgdp_bam_dir}/*${genome_file_id}.srt.bam | bgzip > ${sgdp_vcf_dir}/${genome_file_id}.raw.vcf.gz
#${fermi_program_path}/k8 ${fermi_program_path}/hapdip.js vcfsum -f ${sgdp_vcf_dir}/${genome_file_id}.raw.vcf.gz | bgzip > ${sgdp_vcf_dir}/${genome_file_id}.flt.vcf.gz

#call large variants
#loop over all individuals to call
rm -f ${temp_dir}/all_runs_call_SV_${genome_file_id}.txt

while read -a array
do
    file_id=${array[1]}
    echo -e "bash call_sgdp_SV_variants.sh ${sgdp_bam_dir} ${sgdp_vcf_dir} ${ref_fasta_name} ${file_id} ${genome_file_id}" >> ${temp_dir}/all_runs_call_SV_${genome_file_id}.txt
    #${fermi_program_path}/htsbox abreak -cuf ${ref_file_path} ${sgdp_sam_dir}/${file_name}.sam.gz | gzip -1 > ${sgdp_vcf_dir}/${genome_file_id}.sv.vcf.gz
done<${sgdp_indiv_names_file}


#submit all jobs for calling SV variants
del_dir=/cluster_path/ape_project/deletions_project
totalMem=10G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_call_SV_${genome_file_id}.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_all_call_SV_${genome_file_id}.sh
rm -f ${temp_dir}/submit_all_call_SV_${genome_file_id}.log
qsub -o ${temp_dir}/submit_all_call_SV_${genome_file_id}.log ${temp_dir}/submit_all_call_SV_${genome_file_id}.sh ${temp_dir}/all_runs_call_SV_${genome_file_id}.txt

