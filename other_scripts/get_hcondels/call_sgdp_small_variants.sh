#!/bin/sh
#$ -cwd 
#$ -l h_vmem=30G
#$ -l h_rt=500:00:00
#$ -V
#$ -j y

ref_file_path=$1
sgdp_bam_dir=$2
sgdp_vcf_dir=$3
genome_file_id=$4

fermi_program_path=/cluster_path/bin/fermikit/fermi.kit
${fermi_program_path}/htsbox pileup -cuf ${ref_file_path} ${sgdp_bam_dir}/*${genome_file_id}.srt.bam | bgzip > ${sgdp_vcf_dir}/${genome_file_id}.raw.vcf.gz
${fermi_program_path}/k8 ${fermi_program_path}/hapdip.js vcfsum -f ${sgdp_vcf_dir}/${genome_file_id}.raw.vcf.gz | bgzip > ${sgdp_vcf_dir}/${genome_file_id}.flt.vcf.gz

