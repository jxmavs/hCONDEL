########### merge files across different sequencing lanes first ############

prefix="Rep1-targeted-LOXL2-RNA"
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep1_targeted_loxl2_rna

cat ${data_dir}/${prefix}_S9_L001_R1_001.fastq.gz ${data_dir}/${prefix}_S9_L002_R1_001.fastq.gz ${data_dir}/${prefix}_S9_L003_R1_001.fastq.gz ${data_dir}/${prefix}_S9_L004_R1_001.fastq.gz > ${data_dir}/${prefix}_R1.fastq.gz

prefix="Rep2-targeted-LOXL2-RNA"
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep2_targeted_loxl2_rna

cat ${data_dir}/${prefix}_S10_L001_R1_001.fastq.gz ${data_dir}/${prefix}_S10_L002_R1_001.fastq.gz ${data_dir}/${prefix}_S10_L003_R1_001.fastq.gz ${data_dir}/${prefix}_S10_L004_R1_001.fastq.gz > ${data_dir}/${prefix}_R1.fastq.gz

prefix="Rep1-targeted-LOXL2-RNA"
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep1_targeted_loxl2_rna

cat ${data_dir}/${prefix}_S9_L001_R2_001.fastq.gz ${data_dir}/${prefix}_S9_L002_R2_001.fastq.gz ${data_dir}/${prefix}_S9_L003_R2_001.fastq.gz ${data_dir}/${prefix}_S9_L004_R2_001.fastq.gz > ${data_dir}/${prefix}_R2.fastq.gz

prefix="Rep2-targeted-LOXL2-RNA"
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep2_targeted_loxl2_rna

cat ${data_dir}/${prefix}_S10_L001_R2_001.fastq.gz ${data_dir}/${prefix}_S10_L002_R2_001.fastq.gz ${data_dir}/${prefix}_S10_L003_R2_001.fastq.gz ${data_dir}/${prefix}_S10_L004_R2_001.fastq.gz > ${data_dir}/${prefix}_R2.fastq.gz

########### run GoT ############

run_name=del_1a_2_13_21

script_dir=/cluster_path/ape_project/deletions_project/10X/scripts

run_path=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}
mkdir -p ${run_path}

#download list from here: https://github.com/dan-landau/IronThrone-GoT/tree/master/barcodes10X
whitelist=/cluster_path/bin/IronThrone-GoT-master/barcodes10X/3M-february-2018.txt


fastq_prefix="Rep1-targeted-LOXL2-RNA"
out_prefix=rep1_targeted_loxl2_rna
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/${out_prefix}

config=${run_path}/${run_name}_got_config.txt
fastqR1=${data_dir}/${fastq_prefix}_R1.fastq.gz
fastqR2=${data_dir}/${fastq_prefix}_R2.fastq.gz
outdir=${run_path}/got_output/${out_prefix}
rm -fr ${out_dir}
log=${out_prefix}_log
sample=${out_prefix}

echo -e "${config} ${fastqR1} ${fastqR2} ${whitelist} ${outdir} ${log} ${sample}" > ${run_path}/${run_name}_got_ref.txt


fastq_prefix="Rep2-targeted-LOXL2-RNA"
out_prefix=rep2_targeted_loxl2_rna
data_dir=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/${out_prefix}

config=${run_path}/${run_name}_got_config.txt
fastqR1=${data_dir}/${fastq_prefix}_R1.fastq.gz
fastqR2=${data_dir}/${fastq_prefix}_R2.fastq.gz
outdir=${run_path}/got_output/${out_prefix}
rm -fr ${out_dir}
log=${out_prefix}_log
sample=${out_prefix}

echo -e "${config} ${fastqR1} ${fastqR2} ${whitelist} ${outdir} ${log} ${sample}" >> ${run_path}/${run_name}_got_ref.txt

bash ${script_dir}/GoT_submit.sh ${run_name} ${run_path}/${run_name}_got_ref.txt
