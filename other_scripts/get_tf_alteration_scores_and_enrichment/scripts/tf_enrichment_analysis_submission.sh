####################### run  permutations N times ####################
################## permutation random sampling ################
#${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats is the mismatches
#cons_seq_comp_stat_revised_file=${out_path}/${header_name}_all_stats_hg19_revised_relevant_cols.txt, with hg38 coordinates
#use that to lift to hg38
scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts
roadmap_script_dir=/cluster_path/ape_project/deletions_project/roadmap/scripts
out_path=/cluster_path_temp/tf_enrichment
mkdir -p ${out_path}
wrkdir=${out_path}

tf_enrichment_files_dir=/cluster_path/ape_project/deletions_project/tf_enrichment/files
mkdir -p ${tf_enrichment_files_dir}

#type_file=${tf_enrichment_files_dir}/type_ref_file.txt

#for empirical only
#head -1 ${tf_enrichment_files_dir}/type_ref_file.txt > ${tf_enrichment_files_dir}/type_ref_file.empirical_only.txt
type_file=${tf_enrichment_files_dir}/type_ref_file.empirical_only.txt

species_ref_meta_file=${tf_enrichment_files_dir}/species_ref_input_file.txt
pairwise_comp_file=/cluster_path/ape_project/deletions_project/meme_suite_analysis/pairwise_comp_info.txt

meme_suite_files_path=/cluster_path/ape_project/deletions_project/meme_suite_analysis

cutOffVal=0.0001
tf_db_name=JASPARv2020
memeFileName=${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt

#check that fimo ids are all represented in the meme file
#awk -F"/" '{print $NF}' /cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' | sort -k1,1 > ${out_path}/check_fimo_ids.txt
#join -11 <(sort -k1,1 ${out_path}/check_fimo_ids.txt) -22 <(sort -k2,2 ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt) | wc -l
#rm -f ${out_path}/check_fimo_ids.txt
maxBPWindow=60

fimo_out_path=${out_path}/fimo_output

########## perm parameters ############
maxExtendConsSeqLen=200

lengthBoundPct=0.05
misMatchBoundPct=0.05
GCBoundPct=0.05
logOddsBoundPct=0.05

#for sampling deletions within conserved coordinates, the maximum amount of times to sample
#if we can not sample a deletion site in the conserved coordinate 
sampleNullSeqMaxTrys=500
#the deviation pct from the conserved block of interest if the sequence can not be sampled after maxTrys
sampleNullSeqConsDevPct=0.05

header_name="mostConservedMapped_new_0.001_human_size_filtered"
#perm_header_name="mostCons_all_mismatch_analysis_maxExtendLen_${maxExtendConsSeqLen}_lengthCf_${lengthBoundPct}_mismatchCf_${misMatchBoundPct}_GCCf_${GCBoundPct}_maxTrys_${sampleNullSeqMaxTrys}_devPct_${sampleNullSeqConsDevPct}"
perm_header_name="mostCons_all_mismatch_analysis_maxExtendLen_${maxExtendConsSeqLen}_lengthCf_${lengthBoundPct}_mismatchCf_${misMatchBoundPct}_GCCf_${GCBoundPct}_logOddsCf_${logOddsBoundPct}_maxTrys_${sampleNullSeqMaxTrys}_devPct_${sampleNullSeqConsDevPct}"


#see generate_master_files_seq_tf_enrichment_submission.sh for how this was made
#same files used to generate roadmap permutations
#'seqId', 'chrID','consSeqLen', 'totalMisMatchPct', 'GCPct', 'logOdds'
cons_seq_comp_stat_revised_file=/cluster_path_temp/del_perm_analyses/${header_name}_all_stats_hg38_revised_relevant_cols.txt
orig_seq_comp_stat_revised_file=/cluster_path_temp/del_perm_analyses/MPRA_1_JX_230_BP_36K_oligos_hg38_revised_relevant_cols.txt 

#revise to put hg38 coordinates first
python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin

####### real data parameters ##########

real_data_out_header=hcondel_denylist_filtered

#vim ${tf_enrichment_files_dir}/type_ref_file.txt

########### split meme file ###########

motif_file=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
meme_split_dir=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split
rm -fr ${meme_split_dir}
mkdir -p ${meme_split_dir}
memeSplitFileHeader=${meme_split_dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme
${python_dir}/python ${scriptdir}/seperateMEMEFile.py ${motif_file} ${outFileHeader} ${memeSplitFileHeader}

#generate meme file list
ls -l /cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split | awk -v meme_split_dir=${meme_split_dir} 'NR>1{print meme_split_dir"/"$NF}' > /cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt

################ generate perm data #########################

#run permutations based off of chimp sequence, lift it to hg38? 
outdir=${tf_enrichment_files_dir}/permutations
#rm -fr ${outdir}
mkdir -p ${outdir}

#remember each permutation samples the # in the hcondel set, so if 10032 in hcondel set, then will be 5*10032
numPerm=1000

#${roadmap_script_dir}/gen_null_seq_submit_rand_cons_background_v2_alternate_v2.sh differs from ${roadmap_script_dir}/gen_null_seq_submit_rand_cons_background_v2_alternate.sh in that it accounts for previously sampled sequences
bash ${roadmap_script_dir}/gen_null_seq_submit_rand_cons_background_v2_alternate_v2.sh ${outdir} ${perm_header_name} ${numPerm} ${cons_seq_comp_stat_revised_file} ${orig_seq_comp_stat_revised_file} ${lengthBoundPct} ${misMatchBoundPct} ${GCBoundPct} ${logOddsBoundPct} ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${wrkdir}

head ${wrkdir}/null_seq_submit_${perm_header_name}_all.log
#head /cluster_path_temp/roadmap/null_seq_submit_${perm_header_name}_all.log

#see file output
ls -l ${outdir}/*_list.txt
ls -l ${outdir}/*_all_coord_extended_with_fimo_id.fa | head

outdir=${tf_enrichment_files_dir}/permutations
numPerm=1000
## combine all perm files ##
rm -f ${out_path}/${perm_header_name}_all_coord_extended_with_fimo_id_sorted.fa
for permInd in $(seq 1 ${numPerm})
do
	#revise id
	${python_dir}/python ${scriptdir}/addIdToFasta.py ${outdir}/${permInd}_${perm_header_name}_all_coord_extended_with_fimo_id.fa ${outdir}/${permInd}_${perm_header_name}_all_coord_extended_with_fimo_id.perm_id_added.fa ${permInd}
	cat ${outdir}/${permInd}_${perm_header_name}_all_coord_extended_with_fimo_id.perm_id_added.fa >> ${out_path}/${perm_header_name}_all_coord_extended_with_fimo_id_sorted.fa
done

#double check count first 
#should be 10032*4*1000
#wc -l ${out_path}/${perm_header_name}_all_coord_extended_with_fimo_id.sorted.fa

## split permutations into meta files ##
numCoord=100320
out_split_path=${out_path}/fimo_coord_split
rm -fr ${out_split_path}
mkdir -p ${out_split_path}
${python_dir}/python ${scriptdir}/fimo_split_fasta.py ${out_path}/${perm_header_name}_all_coord_extended_with_fimo_id_sorted.fa ${numCoord} ${out_split_path}/${perm_header_name}_all_coord_extended_with_fimo_id_sorted.fa.

ls -l ${out_split_path}/${perm_header_name}_all_coord_extended_with_fimo_id_sorted.fa.* | awk '{print $NF}' > ${out_path}/${perm_header_name}_all_fimo_file_names.txt


#### double check files - ensure that right number of cooordinates across all files ####
rm -f ${out_path}/perm_files_info.txt
for permInd in $(seq 1 ${numPerm})
do
	wc -l ${outdir}/${permInd}_${perm_header_name}_all_coord_extended_with_fimo_id.fa >> ${out_path}/perm_files_info.txt
done

awk '$1!=40128{print}' ${out_path}/perm_files_info.txt
#40120

while read -a array
do
	fimoFastaIterFile="${array[0]}"
	wc -l ${fimoFastaIterFile}
done<${out_path}/${perm_header_name}_all_fimo_file_names.txt


##################### real data processing #####################

denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt
#no dels that completely overlap conserved block in this set
#awk -F"|" '$9>=$3 && $8<=$2{print}' /cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt | head
awk -F"|" '{print $7"\t"$8"\t"$9"\t+\t"$0}' ${denylist_filtered_ids_file} > ${out_path}/${real_data_out_header}_chimp_coord.txt
awk -F"|" '{print $10"\t"$11"\t"$12"\t"$13"\t"$0}' ${denylist_filtered_ids_file} > ${out_path}/${real_data_out_header}_human_coord.txt


#extend coordinates so that the TF windows can overlap both ways of the site of interest
#extend only coordinates which are common across all species
while read -a array
do
    speciesName="${array[0]}"
    chrSizeFile="${array[2]}"
    #only keep coordinates which can map in all species
    #awk 'NR==FNR{a[$0];next} ($4 in a) {print}' ${out_path}/${real_data_out_header}_commonIDs.txt ${out_path}/${real_data_out_header}_${speciesName}_coord.txt | sort -k4,4 > ${out_path}/${real_data_out_header}_${speciesName}_coord_filtered.txt
    #keep track of inversions
    #python ${scriptdir}/extendCoordinatesCHIPV3b.py ${out_path}/${real_data_out_header}_${speciesName}_coord_filtered.txt ${chrSizeFile} ${maxBPWindow} ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt
	${python_dir}/python ${scriptdir}/extendCoordinatesCHIPV3b.py ${out_path}/${real_data_out_header}_${speciesName}_coord.txt ${chrSizeFile} ${maxBPWindow} ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt
        
    sort -k1,1 -k2,2n ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt > ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_temp.txt
	mv ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_temp.txt ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt


done<${species_ref_meta_file}

rm -f ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id.fa
#get the fasta, append onto master fasta file with all species of interest
while read -a array
do
    speciesName="${array[0]}"
    ref_genome="${array[1]}"
    
    fasta_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/${speciesName}
    
    #create id, add on the species, along with the extended coordinates and the centered coordinate prior to extension
    paste -d"|" <( awk -v ref_genome=${ref_genome} '{print $7"#"ref_genome}' ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt ) <( cut -f1-6 ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt | awk 'BEGIN{OFS="|";} {$1=$1}1'  ) > ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_fimo_id.txt
    
    #revise file header for coordinate
    paste <(cut -f1-3 ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt ) <(cat ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_fimo_id.txt)  > ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.txt
	bedtools getfasta -name -fi ${fasta_ref_dir}/${ref_genome}.fa -bed ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.txt -fo ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa
	
	cat ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa >> ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id.fa


done<${species_ref_meta_file}

#sort combined_fasta_file

awk '!/^>/ { next } { getline seq } { print $0"\t"seq }'  ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id.fa | sort -k1,1 | awk '{print $1"\n"$2}' > ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa

#split fasta file into files so that its numCoord coordinates per file
#this is for computational speedup
#as I have 690 (TFs) * 100 base positions * 10 coordinates * 3 species = 2070000 lines of output

numCoord=1000
#don't delete if not
out_split_path=${out_path}/fimo_coord_split
rm -fr ${out_split_path}
mkdir -p ${out_split_path}
${python_dir}/python ${scriptdir}/fimo_split_fasta.py ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa ${numCoord} ${out_split_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.

ls -l ${out_split_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.* | awk '{print $NF}' > ${out_path}/${real_data_out_header}_all_fimo_file_names.txt

########## submit to generate fimo for all ##############


#reset if running first time
rm -fr ${fimo_out_path}
mkdir -p ${fimo_out_path}

while read -a array
do
	type_header_id=${array[0]}
	seqFileNamesFile=${array[1]}
	memeFileNamesFile=${array[2]}
	
	#run fimo, centered only, only using expressed genes, not checking for encode annotation
	lookAtCenteredOnly=1

	out_header=${type_header_id}_${tf_db_name}

	#remove old files, put 7/9/21 4:18pm
	rm -f ${fimo_out_path}/${out_header}_*
	#v3 takes meme file name as input, same with ${scriptdir}/fimo_run_iter_v4.sh
	#v4 takes an encode ID file to map TFs
	#v5 includes genecard_ID_map.txt to correctly match ids with encode
	bash ${scriptdir}/fimo_run_all_submit_v5_altered.sh ${out_header} ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${memeFileName} ${fimo_out_path} ${seqFileNamesFile} ${memeFileNamesFile} 

	#out_header=tf_enrichment_hcondel_denylist_filtered_random_JASPARv2020
	#head ${fimo_out_path}/${out_header}_fimo_all_runs.txt
	#head ${fimo_out_path}/${out_header}_fimo_all_runs.log
	#ls -l ${fimo_out_path}/*_fimo_all_runs.log

done<${type_file}

#${scriptdir}/fimo_run_iter_v6_altered.sh ${input_file} ${out_header} ${motif_file} ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${memeFileName} ${out_path}

echo -e "${scriptdir}/fimo_run_all_submit_v5_altered.sh ${out_header} ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${memeFileName} ${fimo_out_path} ${seqFileNamesFile} ${memeFileNamesFile}"


##### combine files after finished running ###


#combine files
bash ${scriptdir}/tf_enrichment_combine_submit.sh ${scriptdir} ${out_path} ${pairwise_comp_file} ${tf_enrichment_files_dir} ${type_file} ${tf_db_name} ${fimo_out_path}
#head ${out_path}/tf_enrichment_combine_all_runs.log

#split files
type_header_id=tf_enrichment_hcondel_denylist_filtered_random
randomMemeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
bash ${scriptdir}/tf_enrichment_random_split_submit.sh ${scriptdir} ${out_path} ${tf_db_name} ${fimo_out_path} ${randomMemeFileNamesFile} ${type_header_id} 
#head ${out_path}/tf_enrichment_random_split_all_runs.log 

## check that all tf files have the same # of lines ##
#bash ${scriptdir}/tf_enrichment_perm_file_check.sh 
out_path=/cluster_path/ape_project/deletions_project/tf_enrichment/output
randomMemeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
bash ${scriptdir}/tf_enrichment_random_split_perm_file_check_submit.sh ${scriptdir} ${out_path} ${fimo_out_path} ${randomMemeFileNamesFile}
head ${out_path}/tf_enrichment_random_split_perm_file_check.log

#combine all files
out_path=/cluster_path/ape_project/deletions_project/tf_enrichment/output
randomMemeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
rm -f ${out_path}/all_split_perm_file_list.txt
while read -a array
do
	motif_file=${array[0]}
    tf_name=$( echo ${motif_file} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 2; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' )

	#cat ${fimo_out_path}/perm_final/${tf_name}_file_check.txt >> ${out_path}/all_split_perm_file_list.txt
	cat ${out_path}/${tf_name}_file_check.txt >> ${out_path}/all_split_perm_file_list.txt
done<${randomMemeFileNamesFile}

awk '$1!=10032{print}' ${out_path}/all_split_perm_file_list.txt

### get stats of permutations ###
#old version 
#type_header_id_random=tf_enrichment_hcondel_denylist_filtered_random
#type_header_id_empirical=tf_enrichment_hcondel_denylist_filtered_empirical
#memeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
#bash ${scriptdir}/get_tf_enrichment_perm_stats_old.sh ${type_header_id_random} ${type_header_id_empirical} ${memeFileNamesFile} ${fimo_out_path}

#for random data
out_path=/cluster_path/ape_project/deletions_project/tf_enrichment/output
type_header_id=tf_enrichment_hcondel_denylist_filtered_random
randomMemeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
bash ${scriptdir}/get_tf_enrichment_perm_stats_submit.sh ${scriptdir} ${out_path} ${tf_db_name} ${fimo_out_path} ${randomMemeFileNamesFile} ${type_header_id} 
#head ${out_path}/get_tf_enrichment_perm_stats_all_runs.log

#for empirical data
out_path=/cluster_path/ape_project/deletions_project/tf_enrichment/output
type_header_id=tf_enrichment_hcondel_denylist_filtered_empirical
memeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
bash ${scriptdir}/get_tf_enrichment_empirical_stats.sh ${tf_db_name} ${fimo_out_path} ${memeFileNamesFile} ${type_header_id}

### get stats for each tf/compare with empirical ###
out_path=/cluster_path/ape_project/deletions_project/tf_enrichment/output
type_header_id_empirical=tf_enrichment_hcondel_denylist_filtered_empirical
outFile=${out_path}/tf_enrichment_hcondel_denylist_filtered_final_stats.txt
#perm_summary_stats_file_header=${fimo_out_path}/perm_final/tf_enrichment_hcondel_denylist_filtered_random
perm_summary_stats_file_header=${out_path}/tf_enrichment_hcondel_denylist_filtered_random
memeFileNamesFile=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme_split_files_list.txt
#bash ${scriptdir}/perm_analyze_stats_tf_enrichment.sh ${memeFileNamesFile} ${fimo_out_path}/${type_header_id_empirical}_summary_stats.txt ${outFile} ${perm_summary_stats_file_header}
bash ${scriptdir}/perm_analyze_stats_tf_enrichment.sh ${memeFileNamesFile} ${out_path}/${type_header_id_empirical}_summary_stats.txt ${outFile} ${perm_summary_stats_file_header}

#remove mouse tfs
grep -v -E "Dux|EWSR1-FLI1|Msx3|Rhox11|mix-a" ${out_path}/tf_enrichment_hcondel_denylist_filtered_final_stats.txt > ${out_path}/tf_enrichment_hcondel_denylist_filtered_final_stats.mouse_filtered.txt
	
Dux MA0611.1 DUX DUX
EWSR1-FLI1 MA0149.1 EWSR1-FLI1 EWSR1-FLI1
Msx3 MA0709.1 MSX3 MSX3
Rhox11 MA0629.1 RHOX11 RHOX11
mix-a MA0621.1 MIX-A MIX-A



#remove files, run afterwards
bash ${scriptdir}/tf_enrichment_remove_files_submit.sh ${scriptdir} ${out_path} ${pairwise_comp_file} ${tf_enrichment_files_dir} ${type_file} ${tf_db_name} ${fimo_out_path}

