#4/11/19 nova seq only 
#removed the combining into all lanes from the nova seq files
#and also some comments

#***************************************** novaseq with merging of previous enhancer/barcode data ************************************************************#


############# flash ##############################
#gunzip to temp_dir
#gunzip -c ${file_path}/MPRADel_Enhancer_Barcode_Assoc_1_novaseq_12_15_18_mpradel.fastq.gz > ${temp_dir}/MPRADel_Enhancer_Barcode_Assoc_1_novaseq_12_15_18_mpradel.fastq
#gunzip -c ${file_path}/MPRADel_Enhancer_Barcode_Assoc_2_novaseq_12_15_18_mpradel.fastq.gz > ${temp_dir}/MPRADel_Enhancer_Barcode_Assoc_2_novaseq_12_15_18_mpradel.fastq

#see benchling notes for calculation of fragment size
#but change to average fragment size of library?
#awk 'NR%2==0{print length($0)}' /cluster_path/ape_project/deletions_project/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos.fa | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'
#got 228.828, then add 44 to get 272.828, so 273 bp

read_length=150
fragment_length=273
header_name=enhancer_barcode_1_22_19_merged
temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb
output_dir=${temp_dir}

############### pool all old barcode seqeuencing data, generate ${file_names_file} ###################################

#gzip everything beforehand
data_dir=/cluster_path/ape_project/deletions_project/plasmid_seq/data
header_name=enhancer_barcode_1_22_19_merged
temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb

echo "${data_dir}/hiseq_enhancer_barcode_1_28_17.1.fastq" > ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/hiseq_enhancer_barcode_1_28_17.2.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_1.1.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_1.2.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_2.1.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_2.2.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/nova_seq_test_5_19_18.1.fastq" >> ${data_dir}/${header_name}_gzip_list.txt
echo "${data_dir}/nova_seq_test_5_19_18.2.fastq" >> ${data_dir}/${header_name}_gzip_list.txt


bash gzip_submit.sh ${data_dir}/${header_name}_gzip_list.txt ${header_name} ${temp_dir} 

#run below to generate initial file for the headers for 

file_1_header="Hifi\t${data_dir}/Hifi_1.fastq.gz\t${data_dir}/Hifi_2.fastq.gz"
file_2_header="hiseq_enhancer_barcode_1_28_17\t${data_dir}/hiseq_enhancer_barcode_1_28_17.1.fastq.gz\t${data_dir}/hiseq_enhancer_barcode_1_28_17.2.fastq.gz"
#seperate out lane 1 and lane 2 because lane 2 appeared to have a bit more errors (see Tammy's email on 2/1/18)
file_3_header="hiseq_enhancer_barcode_2_7_18_lane_1\t${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_1.1.fastq.gz\t${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_1.2.fastq.gz"
file_4_header="hiseq_enhancer_barcode_2_7_18_lane_2\t${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_2.1.fastq.gz\t${data_dir}/hiseq_enhancer_barcode_2_7_18_lane_2.2.fastq.gz"
file_5_header="nova_seq_test_5_19_18\t${data_dir}/nova_seq_test_5_19_18.1.fastq.gz\t${data_dir}/nova_seq_test_5_19_18.2.fastq.gz"
file_6_header="MPRADel_Enhancer_Barcode_Assoc_novaseq_12_15_18_mpradel\t${data_dir}/MPRADel_Enhancer_Barcode_Assoc_1_novaseq_12_15_18_mpradel.fastq.gz\t${data_dir}/MPRADel_Enhancer_Barcode_Assoc_2_novaseq_12_15_18_mpradel.fastq.gz"


echo -e "${file_1_header}" > ${data_dir}/${header_name}_combined_file_list.txt
echo -e "${file_2_header}" >> ${data_dir}/${header_name}_combined_file_list.txt
echo -e "${file_3_header}" >> ${data_dir}/${header_name}_combined_file_list.txt
echo -e "${file_4_header}" >> ${data_dir}/${header_name}_combined_file_list.txt
echo -e "${file_5_header}" >> ${data_dir}/${header_name}_combined_file_list.txt
echo -e "${file_6_header}" >> ${data_dir}/${header_name}_combined_file_list.txt

file_names_file=${data_dir}/${header_name}_combined_file_list.txt

#####################################################

#flash_merge_gz.sh allows gzipped as input and output
bash flash_merge_submit.sh ${header_name} ${file_names_file} ${read_length} ${fragment_length} ${output_dir} ${temp_dir} 

#change file names
while read -a array
do
    file_header="${array[0]}"
    mv ${output_dir}/${file_header}.extendedFrags.fastq.gzip ${output_dir}/${file_header}.fastq.gz
done<${file_names_file}
 
#qsub -o ${temp_dir}/${header_name}_flash_merge.txt flash_merge_gz.sh ${read_length} ${fragment_length} ${read_one_file} ${read_two_file} ${output_dir} ${header_name}

#mv ${output_dir}/${header_name}.extendedFrags.fastq.gz ${output_dir}/${header_name}.fastq.gz

#gunzip -c ${output_dir}/${header_name}.extendedFrags.fastq.gz | wc -l
#gunzip -c ${output_dir}/${header_name}.notCombined_1.fastq.gz | wc -l 
#gunzip -c ${output_dir}/${header_name}.notCombined_2.fastq.gz | wc -l 

#qsub -o ${temp_dir}/${header_name}_flash_merge.txt flash_merge.sh ${read_length} ${fragment_length} ${read_one_file} ${read_two_file} ${output_dir} ${header_name}

#mv ${output_dir}/${header_name}.extendedFrags.fastq ${output_dir}/${header_name}.fastq
#gzip ${output_dir}/${header_name}.fastq

################ align files ###################################
#cut off the adapter + adapter with barcode at the ends and align the sequence - I cut simply by truncating the amount of bp
#alternative option is to use say cutadapt
#cutadapt may not be ideal in this case due to the greater error in read 2

header_name=enhancer_barcode_1_22_19_merged
file_names_file=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.txt

out_path=/cluster_path/ape_project/deletions_project/plasmid_seq/analysis
bowtie_ind=/cluster_path/ape_project/deletions_project/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos
#where the curent fastq file is 
file_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb
temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb
#cutAMTBegin=58
#cutAMTEnd=16
#cutAMTEnd removes CGTCAAGCGGCCAGTT (16 bp) adapter sequence from end
#cutAMTBegin removes barcode (20bp barcode sequence) + adapter TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG (38 bp long) 

file_suffix=cutAMTBegin_16_adapterSeqTrimEnd
numCores=5


#v3 outputs algined and unaligned files sepeartely
#it also outputs a trimmed gzippe file and feeds the gzipped file into bowtie2 for alignment
#align all files except MPRADel from novaseq
head -5 ${file_names_file} > /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.sub1.txt
#first 16 adapter bp are automatically trimmed off - hard coded into script, rest of adapter is to be located by string fuzzy search
bash prepare_counts_fastq_pair_merged_submit_v3.sh /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.sub1.txt ${out_path} ${bowtie_ind} ${file_path} ${temp_dir} ${header_name} ${file_suffix} ${numCores}


########## for aligning the MPRADel enhancer/barcode association files, because MPRADel is too big, split MPRADel from nova seq into 1000 smaller files for processing ###################

#takes a list of files and splits it into smaller files
tail -1 ${file_names_file} > /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.sub2.txt
files_to_split_list=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.sub2.txt

#num_lines_split=2500000
num_lines_split=10000000
#file_path is where the original files are located
#split_file_path is where to put the split files
split_out_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb/split
mkdir -p ${split_file_path}
bash split_files_gz_submit.sh ${files_to_split_list} ${file_path} ${split_out_path} ${num_lines_split} ${header_name} ${temp_dir}

#rename split files
#bash split_files_rename.sh ${files_to_split_list} ${split_out_path}
#get new file names file
ls -l ${split_out_path} | awk 'NR>1{print $NF}' | awk -F"/" '{print $NF}' | awk -F"." '{print $1}' > /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_split_file_list.txt
 
split_file_list=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_split_file_list.txt

#submit files that have finished splitting to align
#switch to 2GB and 50 hours
out_path=${temp_dir}
numCores=1
#first 16 adapter bp are automatically trimmed off - hard coded into script, rest of adapter is to be located by string fuzzy search
bash prepare_counts_fastq_pair_merged_submit_v3.sh ${split_file_list} ${out_path} ${bowtie_ind} ${split_out_path} ${temp_dir} ${header_name} ${file_suffix} ${numCores}

#make a new combined file list with the split files?


########### make master file name list with split files #################

#include split files 
grep -v MPRADel_Enhancer_Barcode_Assoc /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.txt | awk '{print $1}' > /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list_with_split_files.txt

cat /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_split_file_list.txt >> /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list_with_split_files.txt


############### get the barcode data #########################################
#min alignment score is -0.6 + -0.6 * L
#http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
#each mismatch is a penalty of -6 for perfect base quality
#so for 230 bp, I can have a min score of -137.4, which corresponds to -137.4/6=22.9 bp mismatches for perfect base qualities 
#mapq filter:  https://biostar.usegalaxy.org/p/18386/
#filter alignment score by having a min of 0.05*230*6=-69?  or -70? this would allow a max of 5% of the 230 bp to have mismatches
#you want to keep all the sequences till this step, because a bad sequence could be associated to a bad barcode, which may be filtered out in prior steps
#so, filter after getting the full barcode table
#setting an alignment score cutoff is probably not a good idea - as shorter sequences would be less penalized

header_name=enhancer_barcode_1_22_19_merged
#cutAMTBegin=58
#cutAMTEnd=16
#associate to barcode
temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb

awk -v file_suffix="cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2" -v file_dir="${temp_dir}" '{print file_dir"/"$1"_"file_suffix".aligned.bam"}' /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list_with_split_files.txt > ${temp_dir}/${header_name}_combined_file_list_bam.aligned.txt

#barcode_fastq_file_name=${temp_dir}/${header_name}.fastq
#get only aligned fastq outputted from bowtie2 to save space
#barcode_fastq_file_name=${temp_dir}/${header_name}_cutAMTBegin_16_adpSeqTrimEnd_trimmed_bowtie2.aligned.fastq.gz
#bam_file_name=${temp_dir}/${header_name}_cutAMTBegin_16_adpSeqTrimEnd_trimmed_bowtie2.aligned.bam
#barcode_len=20

bam_list_file_name=${temp_dir}/${header_name}_combined_file_list_bam.aligned.txt
barcode_assoc_out_file=/cluster_path/ape_project/deletions_project/plasmid_seq/analysis/${header_name}_cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2_barcode_assoc_table.txt
del_assoc_out_file=/cluster_path/ape_project/deletions_project/plasmid_seq/analysis/${header_name}_cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2_del_assoc_table.txt
seq_id_filter_pct=0.95
seq_names_file=/cluster_path/ape_project/deletions_project/barcode_associate/files/MPRA_1_JX_230_BP_36K_oligos_seq_names.txt
 
#barcodeCountAssociationV3.py also directly reads in the bam file from a list of bams to save space and also gets barcode directly from read name and encodes seq names into numbers
#barcodeCountAssociationV4.py is an even more revised version of barcodeCountAssociationV3.py  to save even more memory - all the reads are not stored in one massive list, the statistics are updated as each line is being read 
#convertBarcodeToDelSeqTableV2.py allows filtering barcodes by pct and takes in seq_names_file from encoding done from barcodeCountAssociationV3.py
#dupInfo also stores the total number of other sequences the barcode is linked to, rather than a binary identifier
#NOTE on 1/24/19, I accidentally deleted barcodeCountAssociationV2.py, and replaced with a swap file from last year - it should still be good though
rm -f ${temp_dir}/${header_name}_get_barcode_seq_table.txt
qsub -o ${temp_dir}/${header_name}_get_barcode_seq_table.txt get_barcode_seq_table_v2.sh ${bam_list_file_name} ${barcode_assoc_out_file} ${del_assoc_out_file} ${seq_id_filter_pct} ${seq_names_file}

############# calculate misc stats related to flash/alignment ####################

############# first merge all discarded split files into one file #################

header_name=enhancer_barcode_1_22_19_merged
file_names_file=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.sub2.txt
all_subset_file_names_file=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_split_file_list.txt
suffix="cutAMTBegin_16_adapterSeqTrimEnd_trimmed.discarded.fastq.gz"
split_out_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb
out_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb

temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb

bash cat_merge_submit.sh ${file_names_file} ${all_subset_file_names_file} ${suffix} ${split_out_path} ${out_path} ${header_name} ${temp_dir} 

##################################################################################

header_name=enhancer_barcode_1_22_19_merged
out_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb
split_out_path=/cluster_path_temp/nova_12-15-18/mpradel-ehb/split
all_subset_file_names_file=/cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_split_file_list.txt
temp_dir=/cluster_path_temp/nova_12-15-18/mpradel-ehb
#suffix is for sam file
#cutAMTBegin=58
#cutAMTEnd=16
suffix=cutAMTBegin_16_adapterSeqTrimEnd_trimmed
#get_ehb_stats.sh is original file without any splitting
rm -f ${temp_dir}/get_ehb_stats.txt
qsub -o ${temp_dir}/get_ehb_stats.txt get_ehb_stats_with_split.sh /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.txt ${all_subset_file_names_file} ${out_path} ${header_name} ${suffix}

#bash get_ehb_stats_with_split.sh /cluster_path/ape_project/deletions_project/plasmid_seq/data/${header_name}_combined_file_list.txt ${all_subset_file_names_file} ${split_out_path} ${out_path} ${header_name} ${suffix}

############# stats on barcode table #################
header_name=enhancer_barcode_1_22_19_merged
barcode_assoc_out_file=/cluster_path/ape_project/deletions_project/plasmid_seq/analysis/${header_name}_cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2_barcode_assoc_table.txt

awk '{ if($NF>0){num_dup_barcodes+=$2} else if (($4+$5)/$3 > 0.05) {num_bad_synth_barcodes+=$2} else{a=0}; num_total_reads+=$2 } END { print num_dup_barcodes"\t"num_bad_synth_barcodes"\t"num_total_reads"\t"NR }' ${barcode_assoc_out_file}

#output 2/7/19
#179908101   33793334    2670308421  363446578

#NOTE:
#On 2/21/19:
#I removed MPRADel_Enhancer_Barcode_Assoc fastq, bam and split files because of memory constraints

#the proportion of barcodes bad or dup (8%) amongst all barcodes should be less because you may link more random barcode sequences (causing the denominator may be bigger)
#think venn diagram
############################# 
