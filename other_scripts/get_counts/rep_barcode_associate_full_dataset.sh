############################  includes all previous cell types and NPC data from nextseq 6_21_22 ###########################################################################
del_analysis_path=/cluster_path/ape_project/deletions_project
files_path=/cluster_path_temp/barcode_assoc
mkdir -p ${files_path}

script_path=${del_analysis_path}/barcode_associate/scripts
header_name=all_cell_types_combined_mpradel

#barcode_master_table contains the barcode table from enhancer/tag association
#this is generated from 
barcode_master_table=${del_analysis_path}/plasmid_seq/analysis/enhancer_barcode_1_22_19_merged_cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2_barcode_assoc_table.txt

#where to store temp files
temp_dir=/cluster_path_temp/barcode_assoc
mkdir -p ${temp_dir}

#length of mpra barcode
barcode_length=20

#homology filter for oligo
seq_id_filter_pct=0.95

#oligo_fasta_file is the original fasta file sent to agilent
#the file is also is used to generate seq_names_file
oligo_fasta_file=${del_analysis_path}/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos.fa


########### command to generate seq_names_file #########
seq_names_file=${del_analysis_path}/barcode_associate/files/MPRA_1_JX_230_BP_36K_oligos_seq_names.txt
awk 'NR%2==1{print}' ${oligo_fasta_file} | awk -F">" '{print $2}' > ${seq_names_file}

#file_names_file is generated below
file_names_file=${files_path}/${header_name}_file_names.txt

#deseq_count_file_names contains 3 fields - header name for eventual master count table, file path/name, and condition for deseq
deseq_count_file_names=${files_path}/deseq_${header_name}_file_names.txt

############ make the master condition file/header file names referencing fastq #############

#grab MPRADel files only from the master file list, and keep read 1 only (read 2 is more noisy, look at the fastqc reports)
#${seq_file_info} is from /cluster_path/novaseq_mprau_mpradel/scripts/preprocess_commands.sh

#first column of ${files_path}/${header_name}_file_name_list.txt is the file path
#second column is the output file name
#third column is the lane number
#fourth column is the read number - read 1 or read 2
#fifth and sixth columns are the barcodes

seq_file_info=/cluster_path/novaseq_mprau_mpradel/files/novaseq_12_15_18/novaseq_12_15_18_file_info_local_parsed.txt
grep MPRADel ${seq_file_info} | grep -v MPRADel_Enhancer_Barcode_Assoc | awk '$4==1{print}'  > ${files_path}/${header_name}_file_name_list.txt

#uncomment with nextseq data in the future
cat /cluster_path/nextseq_6_21_22_npc_mpradel/files/nextseq_6_21_22_file_info_local_parsed.txt >> ${files_path}/${header_name}_file_name_list.txt
grep NPC ${files_path}/${header_name}_file_name_list.txt > ${files_path}/${header_name}_file_name_list.npc.txt
 
awk '{print $2}' ${files_path}/${header_name}_file_name_list.txt | sort | uniq > ${files_path}/${header_name}_file_name_list_unique_ids.txt

awk -v header_name=${header_name} '{print $2"_"header_name}' ${files_path}/${header_name}_file_name_list.txt | sort | uniq > ${files_path}/${header_name}_file_name_list_fastq_ids.txt

cat ${files_path}/${header_name}_file_name_list_unique_ids.txt | awk -F"_" '{if(NF==3){print $0"\t"$2} else if(NF==4 || NF==5){print $0"\t"$2"_"$3} else{print $0"\t"$2"_"$3"_"$4"_"$5"_"$6}}' > ${files_path}/${header_name}_cell_type_info.txt

master_condition_file_sorted=${files_path}/${header_name}_cell_type_info.txt

################# combine files from all lanes first ##################################

#gunzip, combine all lanes into one file, then gzip the file
bash ${script_path}/combine_gzip_mult_seq_files_submit.sh ${files_path}/${header_name}_file_name_list.txt ${header_name} ${temp_dir} ${files_path}

#ls -l ${files_path}/*.fastq | awk '{print $NF}' | awk -F"/" '{print $NF}' | awk -F"." '{print $1}' >  ${file_names_file}

############################# generate ${file_names_file} ############################

# command to generate first column of file_names_file - first column is original fastq file, second column is output of count file name #########
#this file may need to be generated manually depeending on preferences
#check if fastq.gz already has header_name, if it has, no need to add it on again
ls -l ${files_path} | grep .fastq.gz | awk '{print $NF}' |  awk -F".gz" '{print $1}' |  awk -v files_path=${files_path} '{print files_path"/"$NF}' > ${files_path}/${header_name}_file_names_col_1.txt
awk -F"/" '{print $NF}' ${files_path}/${header_name}_file_names_col_1.txt | awk -v header_name=${header_name} -v files_path=${files_path} -F"[.]" '{print files_path"/"$1"_seq_count_assoc.txt"}' > ${files_path}/${header_name}_file_names_col_2.txt
paste ${files_path}/${header_name}_file_names_col_1.txt ${files_path}/${header_name}_file_names_col_2.txt > ${files_path}/${header_name}_file_names.txt
rm ${files_path}/${header_name}_file_names_col_1.txt ${files_path}/${header_name}_file_names_col_2.txt


######### main script to get counts for all individual reps ####################

#v2 reads in a gzipped file and also outputs a .stats file for statistics
#it locates the adapter sequence then gets the barcode to account for barcodes less than 20bp
#the read sequence is also reverse complemented
#barcodes which have too many errors in the oligo linked to it are filtered out

bash ${script_path}/get_seq_table_counts_from_barcodes_submit_v2.sh ${barcode_master_table} ${file_names_file} ${seq_names_file} ${temp_dir} ${barcode_length} ${header_name} ${seq_id_filter_pct}

######## after getting counts, merge all count files into one master table table ###############################

#deseq_count_file_names needs to be generated, below are commands of how I did it

header_name=all_cell_types_combined_mpradel
deseq_count_file_names=${temp_dir}/deseq_${header_name}_file_names.txt
readnum=1
rm -f ${deseq_count_file_names}_bar_temp_2
while read -a array
do

    file_ID="${array[0]}"

    echo -e "${files_path}/${file_ID}_${readnum}_${header_name}_seq_count_assoc.txt" >> ${deseq_count_file_names}_bar_temp_2

done<${master_condition_file_sorted}

#first column is id for file, second column is link to actual count file, thrid column contains the deseq condition (i.e. cell type)
paste <(awk '{print $1}' ${master_condition_file_sorted}) <(cat ${deseq_count_file_names}_bar_temp_2 ) > ${deseq_count_file_names}_bar_temp_3
paste <(cat ${deseq_count_file_names}_bar_temp_3) <( awk '{print $2}' ${master_condition_file_sorted} ) > ${deseq_count_file_names}

rm ${deseq_count_file_names}_bar_temp_2 ${deseq_count_file_names}_bar_temp_3


#this generates two files - ${files_path}/${header_name}_deseq.txt, which is the count table per rep - the first column is the sequence name, the second column is the number of unique barcodes, the thrid column is the number of total barcodes
#and ${file_path}/${header_name}_deseq_cond_file.txt, which is the conditions table used for DESEQ
bash ${script_path}/deseq_merge_count_data.sh ${deseq_count_file_names} ${header_name} ${files_path} ${temp_dir}

#get the barcode counts
#this generates ${files_path}/${header_name}_barcode.txt
bash ${script_path}/deseq_merge_count_data_barcodes.sh ${deseq_count_file_names} ${header_name} ${files_path} ${temp_dir}
########### look at barcode analysis stats ################

temp_dir=/cluster_path_temp/barcode_assoc
rm -f ${temp_dir}/${header_name}_barcode_assoc_all_stats.txt
while read -a array
do
    out_file_name="${array[1]}"
    suffix=${header_name}
    file_header=$( echo ${out_file_name} | awk -F"/" '{print $NF}' | awk -F"_"${suffix}  '{print $1}' )
    #get average number of barcodes/reads
    num_avg_uniq_barcodes=$( awk '{total+=$2} END{print total/NR}' ${out_file_name} )
    num_avg_reads=$( awk '{total+=$3} END{print total/NR}' ${out_file_name} )

    #output is:
    #Num_Total_Seq   Num_Total_Barcodes  Num_Barcodes_Not_Linked Num_Bad_Dup_Barcodes_Linked Num_Bad_Synth_Error_Barcodes_Linked Num_Seq_No_Adapter Num_Avg_Barcodes Num_Avg_Reads

    awk -v file_header=${file_header} -v num_avg_uniq_barcodes=${num_avg_uniq_barcodes} -v num_avg_reads=${num_avg_reads} 'NR>1{print file_header"\t"$0"\t"num_avg_uniq_barcodes"\t"num_avg_reads}' ${out_file_name}.stats >>  ${temp_dir}/${header_name}_barcode_assoc_all_stats.txt

done<${file_names_file}

#first column is proportion of good barcodes 
#second column is proportion of barcodes linked
awk '{print ($3-$5-$6)/$3"\t"($3-$4)/$3}' ${temp_dir}/${header_name}_barcode_assoc_all_stats.txt

#save output files
mkdir -p ${del_analysis_path}/barcode_associate/final_files
cp ${temp_dir}/${header_name}_deseq.txt ${del_analysis_path}/barcode_associate/final_files
cp ${temp_dir}/${header_name}_deseq_cond_file.txt ${del_analysis_path}/barcode_associate/final_files
cp ${temp_dir}/${header_name}_barcode_assoc_all_stats.txt ${del_analysis_path}/barcode_associate/final_files
