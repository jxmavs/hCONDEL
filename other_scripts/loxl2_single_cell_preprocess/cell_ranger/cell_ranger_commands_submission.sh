########## create 10X ref file, run cell_ranger ###########

script_dir=/cluster_path/ape_project/deletions_project/10X/scripts

run_name=del_1a_2_13_21
run_path=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}
ref_file=${run_path}/10X_input_${run_name}.txt

id=rep1_rna_out
#download transcriptome from 10X, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
transcriptome=/cluster_path/bin/refdata-gex-GRCh38-2020-A
fastqs=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep1_rna
#can directly link to file off nextseq:
#fastqs=/broad/sabetinextseq/RUNS/210212_NB552060_0107_AH7KKCBGXH/210212_NB552060_0107_AH7KKCBGXH/Alignment_1/20210213_044347/Fastq
sample="Rep1-RNA-1,Rep1-RNA-2,Rep1-RNA-3,Rep1-RNA-4"
expect_cells=7000
localcores=2
localmem=20

echo -e "${id}\t${transcriptome}\t${fastqs}\t${sample}\t${expect_cells}\t${localcores}\t${localmem}" > ${ref_file}

id=rep2_rna_out
#download transcriptome from 10X, https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
transcriptome=/cluster_path/bin/refdata-gex-GRCh38-2020-A
fastqs=/cluster_path/ape_project/deletions_project/10X/runs/${run_name}/files/rep2_rna
#can directly link to file off nextseq:
#fastqs=/broad/sabetinextseq/RUNS/210212_NB552060_0107_AH7KKCBGXH/210212_NB552060_0107_AH7KKCBGXH/Alignment_1/20210213_044347/Fastq
sample="Rep2-RNA-1,Rep2-RNA-2,Rep2-RNA-3,Rep2-RNA-4"
expect_cells=7000
localcores=2
localmem=20

echo -e "${id}\t${transcriptome}\t${fastqs}\t${sample}\t${expect_cells}\t${localcores}\t${localmem}" >> ${ref_file}

bash ${script_dir}/cell_ranger_submit.sh ${run_name} ${run_path} ${ref_file}

###################### aggregate data afterwards ###################### 

bash ${script_dir}/cell_ranger_aggr.sh ${run_name} ${run_path} ${ref_file}
