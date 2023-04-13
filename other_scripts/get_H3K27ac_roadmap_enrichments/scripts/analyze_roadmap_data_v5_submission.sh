################# generate ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_all_stats_subset.txt ##############################

#did I look at conserved block or 200 bp oligo in other - i merged with the big conserved file

fasta_file=/cluster_path/ape_project/deletions_project/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos.fa

roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap

roadmap_script_dir=/cluster_path/ape_project/deletions_project/roadmap/scripts

out_path=/cluster_path_temp/del_perm_analyses
temp_dir=${out_path}

mkdir -p ${out_path}
filterdir=/cluster_path/ape_project/deletions_project/filtering_new

chimp_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp
human_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering

real_data_header=MPRA_1_JX_230_BP_36K_oligos

#added 10/11/21
#denylist filtered ids file, use later for only keeping the ids that passed the denylist filters
denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt

##### additional parameters for permutations sampling random conserved background #####

#this parameter is for script extendCoordinatesOligoDesignConsBackground.py 
maxExtendConsSeqLen=200
#extend coordinates all coordinates to 200 bp length, if over ${maxExtendConsSeqLen} don't extend
#we use this parameter to align the human sequence to the chimp sequence as a control for number of mismatches

perm_rand_cons_background_out_header=panTro4_perm_all_cons_sample_background
lengthBoundPct=0.1
misMatchBoundPct=0.01
GCBoundPct=0.03

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

#### created from gen_master_files_seq_roadmap_enrichment.sh ###

cons_seq_comp_stat_revised_file=${out_path}/${header_name}_all_stats_hg19_revised_relevant_cols.txt
orig_seq_comp_stat_revised_file=${out_path}/${real_data_header}_hg19_revised_relevant_cols.txt 


################### download ######################

 
type=H3K27ac
bash ${roadmap_script_dir}/download_roadmap.sh ${type} ${roadmap_dir}/roadmap_ID_cell_type_list.txt

head ${wrkdir}/download_nonimputed_${type}_all.log
head ${wrkdir}/download_nonimputed_${type}_all.txt
head ${wrkdir}/all_runs_nonimputed_${type}.txt
#download remaining incomplete downloads
#grep -E "E008|E021|E074|E020" ${wrkdir}/all_runs_nonimputed_${type}.txt > ${wrkdir}/all_runs_nonimputed_${type}.remain.txt
grep -E "E003|E017|E026|E037|E045|E056|E061|E068|E078|E080|E104|E109|E124" ${wrkdir}/all_runs_nonimputed_${type}.txt > ${wrkdir}/all_runs_nonimputed_${type}.remain.txt

bash ${wrkdir}/all_runs_nonimputed_${type}.remain.txt

type=H3K27ac
bash ${roadmap_script_dir}/download_roadmap_fc.sh ${type} ${roadmap_dir}/roadmap_ID_cell_type_list.txt

#remove files that are empty(more so for H3K27ac nonimputed)
type=H3K27ac
bash ${roadmap_script_dir}/remove_empty_files.sh ${type}
awk 'FNR==NR { a[$1]; next } ($1 in a)' /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt ${roadmap_dir}/roadmap_ID_cell_type_list.txt > /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/cell_${type}_file_list.txt
paste <(cat /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/cell_${type}_file_list.txt ) <(cut -f2- /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt) > /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_temp.txt
mv /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list_temp.txt /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
#get a list of chromosomes which are in all deletions files



################ decompress ##################################

type=H3K27ac
roadFileType=pval
outdir=/cluster_path_temp/roadmap/nonimputed/${type}
typedir=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt

perm_dir=/cluster_path_temp/roadmap/permutations
deletions_list_file=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_mapped.txt 
mkdir -p ${outdir}
mkdir -p ${perm_dir}
bash ${roadmap_script_dir}/decompressbigwig.sh ${typedir} ${outdir} ${type} ${roadFileType} ${roadmap_ID_file}

wrkdir=/cluster_path_temp/roadmap
head ${wrkdir}/decompress_${roadFileType}_bigwig_${type}_all.log

#decompress fc

type=H3K27ac
roadFileType=fc
outdir=/cluster_path_temp/roadmap/nonimputed/${type}
typedir=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt

perm_dir=/cluster_path_temp/roadmap/permutations
deletions_list_file=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_mapped.txt 
mkdir -p ${outdir}
mkdir -p ${perm_dir}
bash ${roadmap_script_dir}/decompressbigwig.sh ${typedir} ${outdir} ${type} ${roadFileType} ${roadmap_ID_file}

wrkdir=/cluster_path_temp/roadmap
head ${wrkdir}/decompress_${roadFileType}_bigwig_${type}_all.log

#remove huge files to save space
#rm -f broad/hptmp/jxue/roadmap/nonimputed/${type}/*.bedgraph

################## generate stats for real data #########################

del_dir=/cluster_path/ape_project/deletions_project
roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap
wrkdir=${roadmap_dir}/work

#get the relevant cell line info
awk -F"\t" 'NR>3{print $2"\t"$15}' ${roadmap_dir}/jul2013.roadmapData.qc_Consolidated_EpigenomeIDs_summary_Table.tsv | sort -k1,1  > ${roadmap_dir}/roadmap_ID_cell_type_list.txt

#TODO posFile for real data, quick to generate
paste <(awk '{print $NF}' ${out_path}/${real_data_header}_hg19_revised.txt | awk -F'[|]' '{print $10"\t"$11-1"\t"$12}' ) <(awk '{print $NF}'  ${out_path}/${real_data_header}_hg19_revised.txt) | sort -k1,1 -k2,2n > ${out_path}/${real_data_header}_hg19_revised_del_pos.txt
liftOver -minMatch=0.001 ${out_path}/${real_data_header}_hg19_revised_del_pos.txt /cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz ${out_path}/${real_data_header}_hg19_revised_del_pos.hg19_pos_mapped.txt ${out_path}/${real_data_header}_hg19_revised_del_pos.hg19_pos_unmapped.txt

#changed 10/11/21 - change to only keep ids that pass the denylist
join -14 -21 <(sort -k4,4 ${out_path}/${real_data_header}_hg19_revised_del_pos.hg19_pos_mapped.txt) <(sort -k1,1 ${denylist_filtered_ids_file}) | awk '{print $2"\t"$3"\t"$4"\t"$1}' | sort -k1,1 -k2,2n > ${out_path}/${real_data_header}_hg19_revised_del_pos.hg19_pos.txt.temp
mv ${out_path}/${real_data_header}_hg19_revised_del_pos.hg19_pos.txt.temp ${out_path}/${real_data_header}_hg19_revised_pos.txt

posFile=${out_path}/${real_data_header}_hg19_revised_pos.txt

############

type=H3K27ac
outDir=/cluster_path/ape_project/deletions_project/roadmap/work
wigDir=/cluster_path_temp/roadmap/nonimputed/${type}
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
cutOffVal=2
out_header=del_roadmap_nonimputed_rand_cons_background
#posFile=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_primate_filter_mapped_sorted.txt
posFile=${out_path}/${real_data_header}_hg19_revised_pos.txt
bash ${roadmap_script_dir}/extract_roadmap_pval_region_scores_deletion_real_data.sh ${type} ${roadmap_ID_file} ${cutOffVal} ${outDir} ${wigDir} ${out_header} ${posFile}

head ${roadmap_dir}/work/${out_header}_${type}.log 
vim ${wrkdir}/${out_header}_${type}_all_runs.txt
#avg score 
sort -k4,4n ${outDir}/all_${out_header}_${type}_cutOff_${cutOffVal}_stat_out.txt
sort -k5,5n ${outDir}/all_${out_header}_${type}_cutOff_${cutOffVal}_stat_out.txt

#### merge real data scores #####

type=H3K27ac
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
cutOffVal=2
out_header=del_roadmap_nonimputed_rand_cons_background
outDir=/cluster_path/ape_project/deletions_project/roadmap/work
bash ${roadmap_script_dir}/merge_roadmap_score_files_real_data.sh ${type} ${roadmap_ID_file} ${cutOffVal} ${out_header} ${outDir}

head ${outDir}/all_${out_header}_${type}_cutOff_${cutOffVal}_stat_out.txt

#column names are 
#${avgScorePval}\t${propSigCoordPval}\t${totalNumCoordPval}\t${avgScoreFc}\t${propSigCoordFc}\t${totalNumCoordFc}\t${roadmapID}
head ${wrkdir}/all_del_roadmap_nonimputed_rand_cons_background_${type}_cutOff_${cutOffVal}_stat_out.txt
#for looking at cell types which has the most peaks that lie in the region
sort -k2,2n ${wrkdir}/all_del_roadmap_nonimputed_rand_cons_background_${type}_cutOff_${cutOffVal}_stat_out.txt
#for looking at cell types which lie in average highest effect size
sort -k4,4n ${wrkdir}/all_del_roadmap_nonimputed_rand_cons_background_${type}_cutOff_${cutOffVal}_stat_out.txt 



############## remove files after analysis #############
type=H3K27ac
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
cutOffVal=2
out_header=perm_roadmap_rand_cons_background_nonimputed
#out_header=del_roadmap_nonimputed_rand_cons_background
outDir=/cluster_path_temp/roadmap/permutation/files
roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap
#check the output for real data - do I keep it in tmp or in the real directory?
wrkdir=${roadmap_dir}/work

#qsub -o ${out_path}/perm_remove_all_real_data.log perm_remove_all_real_data.sh ${type} ${cutOffVal} ${roadmap_ID_file} ${out_header} ${numPerm} ${wrkdir}

bash ${roadmap_script_dir}/perm_remove_all_real_data.sh ${type} ${cutOffVal} ${roadmap_ID_file} ${out_header} ${wrkdir}

#remove permutation files
bash ${roadmap_script_dir}/perm_remove_all_real_data.sh ${type} ${cutOffVal} ${roadmap_ID_file} ${out_header} ${outDir}
################## permutation random sampling ################

outdir=/cluster_path/ape_project/deletions_project/roadmap/permutations_new
mkdir -p ${outdir}

wrkdir=/cluster_path_temp/del_perm_analyses_temp
mkdir -p ${wrkdir}
numPerm=1000

bash ${roadmap_script_dir}/gen_null_seq_submit_rand_cons_background_v2.sh ${outdir} ${perm_header_name} ${numPerm} ${cons_seq_comp_stat_revised_file} ${orig_seq_comp_stat_revised_file} ${lengthBoundPct} ${misMatchBoundPct} ${GCBoundPct} ${logOddsBoundPct} ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${wrkdir}


head ${wrkdir}/null_seq_submit_${perm_header_name}_all.log
#head /cluster_path_temp/roadmap/null_seq_submit_${out_header}_all.log
echo "bash ${roadmap_script_dir}/gen_null_seq_submit_rand_cons_background_v2.sh ${outdir} ${perm_header_name} ${numPerm} ${cons_seq_comp_stat_revised_file} ${orig_seq_comp_stat_revised_file} ${lengthBoundPct} ${misMatchBoundPct} ${GCBoundPct} ${logOddsBoundPct} ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${wrkdir}"

#see file output
ls -l ${outdir}/*_list.txt

############# generate the statistics for permutations ##################

#type_dir=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}
type=H3K27ac
#roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt

#for testing
#grep "Brain_Dorsolateral_Prefrontal_Cortex" /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt > /cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt.temp
#roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt.temp

roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt

cutOffVal=2
outDir=/cluster_path_temp/roadmap/permutation/files
mkdir -p ${outDir}

wigDir=/cluster_path_temp/roadmap/nonimputed/${type}
#perm_seq_header=hg19_perm_rand_cons_background
perm_seq_header="${perm_header_name}"
out_header=perm_roadmap_rand_cons_background_nonimputed
numPerm=1000
bash ${roadmap_script_dir}/extract_roadmap_pval_region_scores_deletion_perm.sh ${type} ${roadmap_ID_file} ${cutOffVal} ${outDir} ${wigDir} ${perm_seq_header} ${out_header} ${numPerm}

head /cluster_path/ape_project/deletions_project/roadmap/work/${out_header}_${type}_getTypeScoresAll.log

permIter=5
roadmapID=E073
out_header=perm_roadmap_rand_cons_background_nonimputed
head ${outDir}/${out_header}_${permIter}_${roadmapID}_${type}_cutOff_${cutOffVal}_stat_out.txt

############## combine and analyze ################

#combine all permutation files
type=H3K27ac
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
cutOffVal=2
out_header=perm_roadmap_rand_cons_background_nonimputed
numPerm=1000
outDir=/cluster_path_temp/roadmap/permutation/files
roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap
wrkdir=${roadmap_dir}/work

qsub -o ${out_path}/perm_combine_all.log ${roadmap_script_dir}/perm_combine_all.sh ${type} ${cutOffVal} ${roadmap_ID_file} ${out_header} ${numPerm} ${outDir}

del_stat_out_file=${wrkdir}/all_del_roadmap_nonimputed_${type}_cutOff_${cutOffVal}_stat_out.txt

#for testing
#grep "E073" ${wrkdir}/all_del_roadmap_nonimputed_${type}_cutOff_${cutOffVal}_stat_out.txt > ${wrkdir}/all_del_roadmap_nonimputed_${type}_cutOff_${cutOffVal}_stat_out.txt.temp
#del_stat_out_file=${wrkdir}/all_del_roadmap_nonimputed_${type}_cutOff_${cutOffVal}_stat_out.txt.temp

######### perm_analyze_stats_v2 to look at permutation background ##############

del_stat_out_file=${wrkdir}/all_del_roadmap_nonimputed_${type}_cutOff_${cutOffVal}_stat_out.txt

bash ${roadmap_script_dir}/perm_analyze_stats_v2.sh ${type} ${cutOffVal} ${del_stat_out_file} ${roadmap_ID_file} ${out_header} ${numPerm} ${outDir}

#look 
#sort by p-value of proportion significant
sort -t$'\t' -k4,4n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt

#sort by p-value of average score fc (calculated as proportion of perm average fc greater than real average fc )
sort -t$'\t' -k5,5n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt

#sort by real proportion significant
sort -t$'\t' -k8,8n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' 

#sort by real average fc
sort -t$'\t' -k9,9n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' 

#sort by perm proportion significant
sort -t$'\t' -k12,12n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' 

#sort by perm average fc
sort -t$'\t' -k13,13n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' 


#sort by p-value of proportion significant, z-score
awk -F"\t" '{print ($8-$12)/$16}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt  | sort -t$'\t' -k2,2n

awk '{print}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' 

#print proportion significant, avg fold change -real data, as well as permutation and compare
#here, I do see that brain has the most amount of peaks that are significant
#awk -F"[\t]" '{print $1"\t"$2"\t"$9-$8"\t"$13-$12}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | sort -k3,3n
awk -F"[\t]" '{print $1"\t"$2"\t"$8-$12"\t"$9-$13}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | sort -t$'\t' -k3,3n
awk -F"[\t]" '{print $1"\t"$2"\t"$8-$12"\t"$9-$13}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | sort -t$'\t' -k4,4n


#parsed output to look at 
awk -F"[\t]" '{print $1"\t"$2"\t"$4"\t"$8"\t"$12"\t"$16"\t"($8-$12)/$16}' ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt > ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.parsed.look.txt 
mv ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.parsed.look.txt /cluster_path/ape_project/deletions_project/roadmap/processed_files

#look at getDelTypeScoresPerm.sh 
sort -t$'\t' -k8,8n ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats.txt | awk -F"\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$12"\t"$13}' > ${outDir}/all_${out_header}_${type}_cutoff_${cutOffVal}_final_perm_stats_table.txt
${outDir}/all_${out_header}_${roadmapID}_${modType}_cutoff_${cutOffVal}_stats.txt


############## remove files after analysis #############
type=H3K27ac
roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${type}/${type}_file_list.txt
cutOffVal=2
#out_header=perm_roadmap_rand_cons_background_nonimputed
out_header=perm_roadmap_nonimputed
numPerm=1000
outDir=/cluster_path_temp/roadmap/permutation/files
roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap
wrkdir=${roadmap_dir}/work

qsub -o ${out_path}/perm_remove_all.log ${roadmap_script_dir}/perm_remove_all.sh ${type} ${cutOffVal} ${roadmap_ID_file} ${out_header} ${numPerm} ${outDir}
