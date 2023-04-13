#************ set parameters first for all parts ********************#
### initial parameters ###
output_header=magma_analysis_ukbb_all

magma_dir=/cluster_path/ape_project/deletions_project/magma
script_dir=${magma_dir}/scripts

#stores log files and other temporary files
magma_temp_dir=/cluster_path_temp/magma_data
mkdir -p ${magma_temp_dir}

magma_ukbb_dir=${magma_dir}/ukbb

#geneloc_file and additional_meta_file are made in magma_analysis_preprocess_1_window.sh
#four columns containing the gene ID, chromosome, start position and stop position, in that order. It may have a fifth column containing the strand, coded as + for the positive/sense strand and - for the negative/antisense strand.
geneloc_file=${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt

#contains additional variables to condition on for magma - 1st column: gene id, 2nd column: number of conserved blocks, 3rd column: number of conserved bp, 4th column: proportion of gene reg region conserved
additional_meta_file=${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_magma_gene_additional_meta.txt

snp_ref_data=${magma_dir}/ref/g1000_eur_maf_filtered

####### window param ###########
window_size=50000
additional_info="background_all_genes"
settings_header=hCONDEL_window_${window_size}_${additional_info}_hg38_protein_coding

################# remaining files depend on parameter selection above ######################

#${magma_dir}/ref/${settings_header}_gene_overlap_info.txt is made from magma_analysis_preprocess_2_v2_window.r
#first column gene id, second column group
set_file=${magma_dir}/ref/${settings_header}_gene_overlap_info.txt

#********** UKBB preprocess (can skip if already ran once) ***************#
#save UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.csv as tab seperated, unix ending first

awk -F"\t" 'NR>1 {print $1"|"$4"\t"$0}' "${magma_ukbb_dir}/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807_unix_ending.txt" | sort -t$'\t' -k1,1 > ${magma_ukbb_dir}/gwas_v3_with_id.txt

#use GWAS from "primary GWAS" derived from Neale lab 

wget https://www.dropbox.com/s/8vca84rsslgbsua/ukb31063_h2_topline.02Oct2019.tsv.gz?dl=1 -O "${magma_ukbb_dir}/ukb31063_h2_topline.02Oct2019.tsv.gz"

gunzip -c ${magma_ukbb_dir}/ukb31063_h2_topline.02Oct2019.tsv.gz | awk -F"\t" 'NR>1{print $1"|"$28}' | sort -t$'\t' -k1,1 > ${magma_ukbb_dir}/ukbb_single_gwas_to_use.txt

#join to get common GWAS 

join -t $'\t' -11 -21 ${magma_ukbb_dir}/gwas_v3_with_id.txt ${magma_ukbb_dir}/ukbb_single_gwas_to_use.txt | awk -F"\t" '!($6 ~ /v2.tsv.bgz/) {print}' > ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt

#add on additional phenotypic information

wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.both_sexes.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.female.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.male.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.male.tsv.bgz

mv ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.bgz ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.gz
mv ${magma_ukbb_dir}/phenotypes.female.tsv.bgz ${magma_ukbb_dir}/phenotypes.female.tsv.gz
mv ${magma_ukbb_dir}/phenotypes.male.tsv.bgz ${magma_ukbb_dir}/phenotypes.male.tsv.gz
 
gunzip -c ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.gz | awk -F"\t" 'NR>1 {print $1"|both_sexes\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > ${magma_ukbb_dir}/phenotypes.all.with_id.tsv
gunzip -c ${magma_ukbb_dir}/phenotypes.female.tsv.gz | awk -F"\t" 'NR>1 {print $1"|female\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> ${magma_ukbb_dir}/phenotypes.all.with_id.tsv
gunzip -c ${magma_ukbb_dir}/phenotypes.male.tsv.gz | awk -F"\t" 'NR>1 {print $1"|male\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> ${magma_ukbb_dir}/phenotypes.all.with_id.tsv

sort -k1,1 ${magma_ukbb_dir}/phenotypes.all.with_id.tsv > ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv

#fill missing with NA values 
join -eNA -o auto -t $'\t' -a1 -11 -21 ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv | cut -f2- > ${magma_ukbb_dir}/gwas_v3_filtered.final.txt

rm -f ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt

#found out that not all phenotypes are depicted, some are still pending as of 4/10/22
#join -v1 -t $'\t' -11 -21 ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv | head

#****** create template file for GWAS in UKBB and run ********#
rm -f ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.txt

while IFS=$'\t' read -a array
do
	phenotype_code="${array[0]}"
	download_command="${array[5]}"
	download_link=$(echo ${download_command} | awk '{print $2}')
	download_file_header=$(echo ${download_command} | awk '{print $NF}' | awk -F"." '{for (i = 1; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}')
	
	echo -e "bash ${magma_dir}/scripts/magma_analyze_ukbb.sh\t${phenotype_code}\t${download_link}\t${download_file_header}\t${magma_temp_dir}\t${settings_header}\t${geneloc_file}\t${snp_ref_data}\t${set_file}\t${additional_meta_file}" >> ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.txt

done<${magma_ukbb_dir}/gwas_v3_filtered.final.txt

############# magma results ################

totalMem=5G
totalNumTasks=$(wc -l ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${script_dir}/runTasksEval.sh > ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.sh

rm -f ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.log
qsub -o ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.log ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.sh ${magma_temp_dir}/${settings_header}_get_magma_results_ukbb.txt


#********** combine magma results after running ***************#

####### create meta file first ###########

echo -e "phenotype_code\tphenotype\tsex" > ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.txt
awk -F"\t" '{print $1"\t"$2"\t"$4}' ${magma_ukbb_dir}/gwas_v3_filtered.final.txt >> ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.txt

echo -e "n_non_missing\tn_missing\tn_controls\tn_cases" > ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.2.txt
awk -F"\t" '{print $12"\t"$13"\t"$14"\t"$15}' ${magma_ukbb_dir}/gwas_v3_filtered.final.txt >> ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.2.txt

#get N, do this for now because some phenotype samples don't have N
rm -f ${magma_dir}/${output_header}_${settings_header}_N_all.txt
echo "n_complete_samples" > ${magma_dir}/${output_header}_${settings_header}_N_all.txt

while read -a array
do
	phenotype_code="${array[0]}"
	
	#sanity check that there are not multiple N's in GWAS file
	num_N=$(wc -l ${magma_temp_dir}/${phenotype_code}/${phenotype_code}_magma_analysis_${settings_header}_N.txt | awk '{print $1}')
	
	if [[ ${num_N} -gt 1 ]]
	then
		echo "multiple N"
	fi
	
	awk '{print $1}' ${magma_temp_dir}/${phenotype_code}/${phenotype_code}_magma_analysis_${settings_header}_N.txt >> ${magma_dir}/${output_header}_${settings_header}_N_all.txt
	
done<${magma_ukbb_dir}/gwas_v3_filtered.final.txt

paste <(cat ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.txt) <(cat ${magma_dir}/${output_header}_${settings_header}_N_all.txt) <(cat ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.2.txt) > ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.txt

rm -f ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.txt
rm -f ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.temp.2.txt

############## pool files ################################

rm -f ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

echo -e "beta\tbeta_std\tse\tpval" > ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

while read -a array
do
	phenotype_code="${array[0]}"
	grep -v "#" ${magma_temp_dir}/${phenotype_code}/${phenotype_code}_magma_analysis_${settings_header}_gs_condition_cons_prop.gsa.out | awk '$1=="Overlap" && $3=="2"{print $5"\t"$6"\t"$7"\t"$8}' >> ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt

done<${magma_ukbb_dir}/gwas_v3_filtered.final.txt

#add on phenotype code and phenotype name
paste <(cat ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.txt) <(cat ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_stat_out.txt) > ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_summary.txt

sort -t$'\t' -k12,12g ${magma_dir}/${output_header}_${settings_header}_condition_cons_prop_summary.txt | head 

######## MISC ###########

#add on additional phenotypic information

wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.both_sexes.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.female.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.female.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.male.tsv.bgz -O ${magma_ukbb_dir}/phenotypes.male.tsv.bgz

mv ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.bgz ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.gz
mv ${magma_ukbb_dir}/phenotypes.female.tsv.bgz ${magma_ukbb_dir}/phenotypes.female.tsv.gz
mv ${magma_ukbb_dir}/phenotypes.male.tsv.bgz ${magma_ukbb_dir}/phenotypes.male.tsv.gz
 
gunzip -c ${magma_ukbb_dir}/phenotypes.both_sexes.tsv.gz | awk -F"\t" 'NR>1 {print $1"|both_sexes\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > ${magma_ukbb_dir}/phenotypes.all.with_id.tsv
gunzip -c ${magma_ukbb_dir}/phenotypes.female.tsv.gz | awk -F"\t" 'NR>1 {print $1"|female\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> ${magma_ukbb_dir}/phenotypes.all.with_id.tsv
gunzip -c ${magma_ukbb_dir}/phenotypes.male.tsv.gz | awk -F"\t" 'NR>1 {print $1"|male\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> ${magma_ukbb_dir}/phenotypes.all.with_id.tsv

sort -k1,1 ${magma_ukbb_dir}/phenotypes.all.with_id.tsv > ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv

join -t $'\t' -11 -21 ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv | cut -f2- > ${magma_ukbb_dir}/gwas_v3_filtered.final.txt

#found out that not all phenotypes are depicted, some are still pending as of 4/10/22
join -v1 -t $'\t' -11 -21 ${magma_ukbb_dir}/gwas_v3_filtered.temp.txt ${magma_ukbb_dir}/phenotypes.all.with_id.sorted.tsv | head

####### create meta file for later ###########

echo -e "phenotype_code\tphenotype\tN" > ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.txt
awk '{print $1"\t"$2"\t"$9}' ${magma_ukbb_dir}/gwas_v3_filtered.final.txt > ${magma_ukbb_dir}/ukbb_magma_summary_df_initial_info.txt


