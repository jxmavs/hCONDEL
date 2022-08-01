
sgdp_bam_dir=$1
sgdp_vcf_dir=$2
ref_file_path=$3
file_id=$4
genome_file_id=$5

fermi_program_path=/cluster_path/bin/fermikit/fermi.kit

#call SV

#${fermi_program_path}/htsbox abreak -bcuf ${ref_file_path} ${sgdp_bam_dir}/${file_id}_${genome_file_id}.srt.bam | gzip -1 > ${sgdp_vcf_dir}/${file_id}_${genome_file_id}.sv.vcf.gz
#s flag is for alignment score 
#q flag is for mapQ quality
${fermi_program_path}/htsbox abreak -c -s0 -b -q0 -f ${ref_file_path} ${sgdp_bam_dir}/${file_id}_${genome_file_id}.srt.bam | gzip -1 > ${sgdp_vcf_dir}/${file_id}_${genome_file_id}.sv.vcf.gz
echo "${fermi_program_path}/htsbox abreak -c -s0 -b -q0 -f ${ref_file_path} ${sgdp_bam_dir}/${file_id}_${genome_file_id}.srt.bam | gzip -1 > ${sgdp_vcf_dir}/${file_id}_${genome_file_id}.sv.vcf.gz"

#${fermi_program_path}/htsbox abreak -c -s0 -b -q40 -f ${ref_file_path} ${sgdp_bam_dir}/${file_id}_${genome_file_id}.srt.bam | gzip -1 > ${sgdp_vcf_dir}/${file_id}_${genome_file_id}.sv.vcf.gz
