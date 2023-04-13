##### generate meta files for roadmap enrichment, which is used in analyze_roadmap_data_v5_submission.sh #####

fasta_file=/cluster_path/ape_project/deletions_project/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos.fa

out_path=/cluster_path_temp/del_perm_analyses

real_data_header=MPRA_1_JX_230_BP_36K_oligos

#added 10/11/21
#denylist filtered ids file, use later for only keeping the ids that passed the denylist filters
denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt

header_name="mostConservedMapped_new_0.001_human_size_filtered"

#${cons_seq_comp_stat_file} and ${orig_seq_comp_stat_file} are pre-generated files from gen_mismatch_conservation_stats.sh

#TODO: change file names below to be more readable?

cons_seq_comp_stat_file=${out_path}/${header_name}_all_stats.txt
orig_seq_comp_stat_file=${out_path}/${real_data_header}.cons_stats.txt

#************ for generating the meta file for randomly sampling conserved blocks, matching based on cons. block size ************#

# also match on number of differences with human, as well as number of deletion sitting in each block #

#liftOver all the conserved sequences to hg19 instead
#$4"\t"$5"\t"$6 are hg38 coordinates
awk '{print $1}' ${cons_seq_comp_stat_file} | awk -F'[|]' '{print $4"\t"$5"\t"$6"\t"$0}' > ${out_path}/${header_name}_all_stats_hg38_coord.txt
liftOver -minMatch=0.001 ${out_path}/${header_name}_all_stats_hg38_coord.txt /cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz ${out_path}/${header_name}_all_stats_hg38_coord.mapped.txt ${out_path}/${header_name}_all_stats_hg38_coord.unmapped.txt

#get the statistics back after lifting over
join -t$'\t' -14 ${out_path}/${header_name}_all_stats_hg38_coord.mapped.txt -21 ${cons_seq_comp_stat_file} > ${out_path}/${header_name}_all_stats_hg19.txt

#revise the file now to replace hg38 coord with hg19 in the id
#first 3 coordinates are panTro4 cons coordinates, then next 3 are hg38 cons coordinates, followed by columns 7,8,9 being hg19 coordinates
#output oldid followed by new id
paste <( awk -F'[|\t]' '{print $1"|"$2"|"$3"|"$4"|"$5"|"$6"\t"$1"|"$2"|"$3"|"$7"|"$8"|"$9}' ${out_path}/${header_name}_all_stats_hg19.txt ) <( cut -f5- ${out_path}/${header_name}_all_stats_hg19.txt ) > ${out_path}/${header_name}_all_stats_hg19_temp.txt


### join to get the statistics for the original deletions ###
paste <(awk '{print $1}' ${out_path}/${real_data_header}_id.txt | awk -F'[|]' '{print $1"|"$2"|"$3"|"$4"|"$5"|"$6}') <(awk '{print $1}' ${out_path}/${real_data_header}_id.txt)  >  ${out_path}/${real_data_header}_temp_ids.txt 

#for the first 3 fields, the first id is the original conserved chimp/human coordinates, the next is the full ID, the final is the new id with hg19 
join -t$'\t' -11 ${out_path}/${real_data_header}_temp_ids.txt -21 ${out_path}/${header_name}_all_stats_hg19_temp.txt | cut -f2-  > ${out_path}/${real_data_header}_hg19_revised_temp.txt 

#put the original hcondel ID also at the end for use later
paste <(cut -f2- ${out_path}/${real_data_header}_hg19_revised_temp.txt ) <( awk '{print $1}' ${out_path}/${real_data_header}_hg19_revised_temp.txt  ) > ${out_path}/${real_data_header}_hg19_revised.txt 

#the revised file contains the human hg38 coordinates replaced with the human hg19 coordinates in the id as the first column
cut -f2- ${out_path}/${header_name}_all_stats_hg19_temp.txt > ${out_path}/${header_name}_all_stats_hg19_revised.txt

rm -f ${out_path}/${header_name}_all_stats_hg19_temp.txt 
rm -f ${out_path}/${real_data_header}_temp_ids.txt

#keep only relevant columns for python run
cons_seq_comp_stat_revised_file=${out_path}/${header_name}_all_stats_hg19_revised_relevant_cols.txt
orig_seq_comp_stat_revised_file=${out_path}/${real_data_header}_hg19_revised_relevant_cols.txt 
 
#output 
#str(totalNumMismatch)+"\t"+str(totalNumGapOpen)+"\t"+str(totalNumGaps)+"\t"+str(numUnalignedBp)+"\t"+str(numUnalignedBp+totalNumMismatch+totalNumGaps)
#'seqId', 'chrID','consSeqLen', 'totalMisMatchPct', 'GCPct', 'logOdds'
awk '{print $1"\t"$2"\t"$3"\t"$8/$3"\t"$11"\t"$17}' ${out_path}/${header_name}_all_stats_hg19_revised.txt > ${cons_seq_comp_stat_revised_file}
awk '{print $1"\t"$2"\t"$3"\t"$8/$3"\t"$11"\t"$17"\t"$NF}' ${out_path}/${real_data_header}_hg19_revised.txt > ${orig_seq_comp_stat_revised_file}

#changed 10/11/21 - change to only keep ids that pass the denylist
awk '{print $1"\t"$2"\t"$3"\t"$8/$3"\t"$11"\t"$17"\t"$NF}' ${out_path}/${real_data_header}_hg19_revised.txt > ${orig_seq_comp_stat_revised_file}.temp
join -17 -21 <(sort -k7,7 ${orig_seq_comp_stat_revised_file}.temp) <(sort -k1,1 ${denylist_filtered_ids_file}) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$1}'  >  ${orig_seq_comp_stat_revised_file}

