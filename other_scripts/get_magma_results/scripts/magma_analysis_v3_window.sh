#note in v3, I separated the preprocess commands to magma_analysis_prerprocess_1.sh

#*********************** initial parameters ***********************#
output_header=magma_analysis

magma_dir=/cluster_path/ape_project/deletions_project/magma
script_dir=${magma_dir}/scripts

#stores log files and other temporary files
magma_temp_dir=/cluster_path/ape_project/deletions_project/magma/temp
mkdir -p ${magma_temp_dir}

#four columns containing the gene ID, chromosome, start position and stop position, in that order. It may have a fifth column containing the strand, coded as + for the positive/sense strand and - for the negative/antisense strand.
geneloc_file=${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt

#contains additional variables to condition on for magma - 1st column: gene id, 2nd column: number of conserved blocks, 3rd column: number of conserved bp, 4th column: proportion of gene reg region conserved
additional_meta_file=${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_magma_gene_additional_meta.txt

snp_ref_data=${magma_dir}/ref/g1000_eur_maf_filtered

#*********************** set initial info running - pick one of these parameter settings ***********************#

####### window param ###########
window_size=50000
additional_info="background_all_genes"
settings_header=hCONDEL_window_${window_size}_${additional_info}_hg38_protein_coding

####### window param random ###########
window_size=50000
additional_info="background_all_genes"
settings_header=hCONDEL_window_${window_size}_${additional_info}_random_matched_hg38_protein_coding

################# remaining files depend on parameter selection above ######################

#${magma_dir}/ref/${settings_header}_gene_overlap_info.txt is made from magma_analysis_preprocess_2_v2.r
#first column gene id, second column group
set_file=${magma_dir}/ref/${settings_header}_gene_overlap_info.txt

rm -f ${magma_temp_dir}/${settings_header}_get_magma_results.txt

####################################### Schizophrenia ##########################################
#https://pubmed.ncbi.nlm.nih.gov/29906448/
#https://www.cell.com/cell/fulltext/S0092-8674(18)30658-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418306585%3Fshowall%3Dtrue#secsectitle0180
#used EU panel for LD for Regional SNP-heritability estimation
#all variants align to human genome build 19 (hg19)
#file location
#https://figshare.com/articles/dataset/cdg2018-bip-scz/14672019
#33,426 SCZ cases and 32,541 SCZ controls, added is 65,967

analysis_header=BDSWGPGC_2018_magma_analysis
work_dir=${magma_dir}/schiz
mkdir -p ${work_dir}
snploc_file_initial=${work_dir}/sczvscont-sumstat.gz

#generate snploc file
#gunzip -c ${snploc_file_initial} | awk 'NR>1{print $1"_"$3"_"$4"_"$5"\t"$1"\t"$3}' > ${work_dir}/${analysis_header}_magma_snploc_input.txt
#generate pval_file
#gunzip -c ${snploc_file_initial} | awk 'NR>1{print $1"_"$3"_"$4"_"$5"\t"$11}' > ${work_dir}/${analysis_header}_magma_pval_input.txt

#a header is not allowed, file must have three columns containing the SNP ID, chromosome and base pair position
snploc_file=${work_dir}/${analysis_header}_magma_snploc_input.txt

annot_prefix=${work_dir}/${analysis_header}_${settings_header}_annot

#file with previously computed SNP p-values (in columns ÔSNPÕ and ÔPÕ in the file [PVAL_FILE])
pval_file=${work_dir}/${analysis_header}_magma_pval_input.txt
gene_prefix=${work_dir}/${analysis_header}_${settings_header}_gene
N=65967

gs_prefix=${work_dir}/${analysis_header}_${settings_header}_gs

#remove output files from potential previous runs
rm -f ${gene_prefix}*
rm -f ${gs_prefix}*

echo -e "bash /cluster_path/ape_project/deletions_project/magma/scripts/magma_analyze.sh\t${snploc_file}\t${geneloc_file}\t${annot_prefix}\t${snp_ref_data}\t${pval_file}\t${gene_prefix}\t${N}\t${set_file}\t${additional_meta_file}\t${gs_prefix}" >> ${magma_temp_dir}/${settings_header}_get_magma_results.txt
#bash ${magma_temp_dir}/${analysis_header}_get_magma_results.txt

####################################### Intelligence ##########################################
#https://www.nature.com/articles/s41588-018-0152-6
#https://ctg.cncr.nl/software/summary_statistics
#wget https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip -P ${magma_dir}/intel
#unzip ${snploc_file_initial}
#269,867
#hg19

analysis_header=savage_2018_magma_analysis
work_dir=${magma_dir}/intel
mkdir -p ${work_dir}
snploc_file_initial=${work_dir}/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt

#generate snploc file
#cat ${snploc_file_initial} | awk 'NR>1{print $3"_"$4"_"toupper($5)"_"toupper($6)"\t"$3"\t"$4}' > ${work_dir}/${analysis_header}_magma_snploc_input.txt
#generate pval_file
#cat ${snploc_file_initial} | awk 'NR>1{print $3"_"$4"_"toupper($5)"_"toupper($6)"\t"$11"\t"$12}' > ${work_dir}/${analysis_header}_magma_pval_input.txt

#a header is not allowed, file must have three columns containing the SNP ID, chromosome and base pair position
snploc_file=${work_dir}/${analysis_header}_magma_snploc_input.txt

annot_prefix=${work_dir}/${analysis_header}_${settings_header}_annot

#file with previously computed SNP p-values (in columns ÔSNPÕ and ÔPÕ in the file [PVAL_FILE])
pval_file=${work_dir}/${analysis_header}_magma_pval_input.txt
gene_prefix=${work_dir}/${analysis_header}_${settings_header}_gene
#N=269867
N="file"

gs_prefix=${work_dir}/${analysis_header}_${settings_header}_gs

#remove output files from potential previous runs
rm -f ${gene_prefix}*
rm -f ${gs_prefix}*

echo -e "bash /cluster_path/ape_project/deletions_project/magma/scripts/magma_analyze.sh\t${snploc_file}\t${geneloc_file}\t${annot_prefix}\t${snp_ref_data}\t${pval_file}\t${gene_prefix}\t${N}\t${set_file}\t${additional_meta_file}\t${gs_prefix}" >> ${magma_temp_dir}/${settings_header}_get_magma_results.txt
#bash ${magma_temp_dir}/${analysis_header}_get_magma_results.txt

####################################### Depression 2 ##########################################
#https://www.nature.com/articles/s41588-018-0090-3
#https://figshare.com/articles/dataset/mdd2018/14672085

#tested 2 depression datasets

#first number cases, second number controls
#total
#135,458+344,901=480359
#23andMe
#75,607+231,747=307354
#480359-307354=173005

analysis_header=wray_2018_magma_analysis
work_dir=${magma_dir}/dep2
mkdir -p ${work_dir}
snploc_file_initial=${work_dir}/MDD2018_ex23andMe.gz

#generate snploc file
#gunzip -c ${snploc_file_initial} | awk 'NR>1{print $1"_"$3"_"$4"_"$5"\t"$1"\t"$3}' > ${work_dir}/${analysis_header}_magma_snploc_input.txt
#generate pval_file
#gunzip -c ${snploc_file_initial} | awk 'NR>1{print $1"_"$3"_"$4"_"$5"\t"$11"\t"$17+$18}' > ${work_dir}/${analysis_header}_magma_pval_input.txt

#a header is not allowed, file must have three columns containing the SNP ID, chromosome and base pair position
snploc_file=${work_dir}/${analysis_header}_magma_snploc_input.txt

annot_prefix=${work_dir}/${analysis_header}_${settings_header}_annot

#file with previously computed SNP p-values (in columns ÔSNPÕ and ÔPÕ in the file [PVAL_FILE])
pval_file=${work_dir}/${analysis_header}_magma_pval_input.txt
gene_prefix=${work_dir}/${analysis_header}_${settings_header}_gene
N="file"

gs_prefix=${work_dir}/${analysis_header}_${settings_header}_gs

#remove output files from potential previous runs
rm -f ${gene_prefix}*
rm -f ${gs_prefix}*

echo -e "bash /cluster_path/ape_project/deletions_project/magma/scripts/magma_analyze.sh\t${snploc_file}\t${geneloc_file}\t${annot_prefix}\t${snp_ref_data}\t${pval_file}\t${gene_prefix}\t${N}\t${set_file}\t${additional_meta_file}\t${gs_prefix}" >> ${magma_temp_dir}/${settings_header}_get_magma_results.txt
#bash ${magma_temp_dir}/${analysis_header}_get_magma_results.txt

####################################### Bipolar ##########################################

#https://www.nature.com/articles/s41588-021-00857-4#data-availability
#https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594

analysis_header=mullins_2021_magma_analysis
work_dir=${magma_dir}/bipolar
mkdir -p ${work_dir}
snploc_file_initial=${work_dir}/pgc-bip2021-all.vcf.tsv.gz

#generate snploc file
#gunzip -c ${snploc_file_initial} | grep -v "#" | awk '{print $1"_"$2"_"$4"_"$5"\t"$1"\t"$2}' > ${work_dir}/${analysis_header}_magma_snploc_input.txt
#generate pval_file
#gunzip -c ${snploc_file_initial} | grep -v "#" | awk '{print $1"_"$2"_"$4"_"$5"\t"$8"\t"$13*2}' > ${work_dir}/${analysis_header}_magma_pval_input.txt
#41,917 individuals with BD (cases) and 371,549 controls of European descent

#a header is not allowed, file must have three columns containing the SNP ID, chromosome and base pair position
snploc_file=${work_dir}/${analysis_header}_magma_snploc_input.txt

annot_prefix=${work_dir}/${analysis_header}_${settings_header}_annot

#file with previously computed SNP p-values (in columns ÔSNPÕ and ÔPÕ in the file [PVAL_FILE])
pval_file=${work_dir}/${analysis_header}_magma_pval_input.txt
gene_prefix=${work_dir}/${analysis_header}_${settings_header}_gene
#N=413466
N="file"

gs_prefix=${work_dir}/${analysis_header}_${settings_header}_gs

#remove output files from potential previous runs
rm -f ${gene_prefix}*
rm -f ${gs_prefix}*

echo -e "bash /cluster_path/ape_project/deletions_project/magma/scripts/magma_analyze.sh\t${snploc_file}\t${geneloc_file}\t${annot_prefix}\t${snp_ref_data}\t${pval_file}\t${gene_prefix}\t${N}\t${set_file}\t${additional_meta_file}\t${gs_prefix}" >> ${magma_temp_dir}/${settings_header}_get_magma_results.txt
#bash ${magma_temp_dir}/${analysis_header}_get_magma_results.txt

############# run all magma commands ################

totalMem=5G
totalNumTasks=$(wc -l ${magma_temp_dir}/${settings_header}_get_magma_results.txt  | awk '{print $1}')
totalTime=20:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${script_dir}/runTasksEval.sh > ${magma_temp_dir}/${settings_header}_get_magma_results.sh

rm -f ${magma_temp_dir}/${settings_header}_get_magma_results.log
qsub -o ${magma_temp_dir}/${settings_header}_get_magma_results.log ${magma_temp_dir}/${settings_header}_get_magma_results.sh ${magma_temp_dir}/${settings_header}_get_magma_results.txt
#head ${magma_temp_dir}/${settings_header}_get_magma_results.log

#*********************** combine output files after running ***********************#
#run each of the settings below one at a time to get the file for each setting

####### window param ###########
window_size=50000
additional_info="background_all_genes"
settings_header=hCONDEL_window_${window_size}_${additional_info}_hg38_protein_coding

####### window param random ###########
window_size=50000
additional_info="background_all_genes"
settings_header=hCONDEL_window_${window_size}_${additional_info}_random_matched_hg38_protein_coding

####### create output stat file from analyses #################
rm -f ${magma_dir}/${output_header}_${settings_header}_stat_out.txt
rm -f ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

echo -e "beta\tbeta_std\tse\tpval" > ${magma_dir}/${output_header}_${settings_header}_stat_out.txt
echo -e "beta\tbeta_std\tse\tpval" > ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

while read -a array
do
	file_header="${array[0]}"
	grep -v "#" ${magma_dir}/${file_header}_magma_analysis_${settings_header}_gs.gsa.out | awk '$1=="Overlap" {print $4"\t"$5"\t"$6"\t"$7}' >> ${magma_dir}/${output_header}_${settings_header}_stat_out.txt
	grep -v "#" ${magma_dir}/${file_header}_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out | awk '$1=="Overlap" && $3=="2"{print $5"\t"$6"\t"$7"\t"$8}' >> ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

done<${magma_dir}/magma_analysis_file_headers_window.txt

paste <(cat ${magma_dir}/magma_summary_df_initial_info_window.txt) <(cat ${magma_dir}/${output_header}_${settings_header}_stat_out.txt) > ${magma_dir}/${output_header}_${settings_header}_summary.txt
paste <(cat ${magma_dir}/magma_summary_df_initial_info_window.txt) <(cat ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt) > ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_summary.txt


#*********************** misc ***********************#

head -20 ${magma_dir}/schiz/BDSWGPGC_2018_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out

head -20 ${magma_dir}/intel/savage_2018_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out

head -20 ${magma_dir}/dep2/wray_2018_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out

head -20 ${magma_dir}/bipolar/mullins_2021_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out
