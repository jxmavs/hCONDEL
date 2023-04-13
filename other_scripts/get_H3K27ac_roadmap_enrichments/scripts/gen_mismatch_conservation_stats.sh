################# generate ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_all_stats_subset.txt ##############################
#use BLAST to get stats
use BLAST+

#did I look at conserved block or 200 bp oligo in other - i merged with the big conserved file

fasta_file=/cluster_path/ape_project/deletions_project/mpra_1_set/MPRA_1_JX_230_BP_36K_oligos.fa

script_dir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts
script_dir_2=/cluster_path/ape_project/deletions_project/LCL_chip/scripts

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


##########################################################################


awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt > ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_tab_sep.txt

#extend coordinates all coordinates to 200 bp length, if over ${maxExtendConsSeqLen} don't extend
#first 6 coordinates contains the new extended sequence after adding to ensure that the min of chimp seq/human seq is 200 bp

python ${script_dir_2}/extendCoordinatesOligoDesignConsBackground.py ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_tab_sep.txt ${chimp_ref_dir}/panTro4.chrom.sizes ${human_ref_dir}/hg38.chrom.sizes ${maxExtendConsSeqLen} ${out_path}/${header_name}_extended_bp.txt

inputFile=${out_path}/${header_name}_extended_bp.txt

#the id will be chimp coord, human coord
paste -d "\t" <(awk '{print $1"\t"$2"\t"$3}' ${inputFile} ) <( cut -f7- ${inputFile} | awk 'BEGIN{OFS="|";} {$1=$1}1') | sort -k1,1 -k2,2n | uniq > ${out_path}/${header_name}_chimp_coord.txt 

paste -d "\t" <(awk '{print $4"\t"$5"\t"$6}' ${inputFile} ) <( cut -f7- ${inputFile} | awk 'BEGIN{OFS="|";} {$1=$1}1' ) | sort -k1,1 -k2,2n | uniq > ${out_path}/${header_name}_human_coord.txt 

#getting fasta sequences for creating the oligos - instead of reading whole ref genome, it serves as a major computational speedup
bedtools getfasta -name -fi ${chimp_ref_dir}/panTro4.fa -bed ${out_path}/${header_name}_chimp_coord.txt -fo ${out_path}/${header_name}_chimp_coord_panTro4.fa
bedtools getfasta -name -fi ${human_ref_dir}/hg38.fa -bed ${out_path}/${header_name}_human_coord.txt -fo ${out_path}/${header_name}_human_coord_hg38.fa

#add on id the fasta sequences
chimp_fasta_file=${out_path}/${header_name}_chimp_coord_panTro4.fa
awk '!/^>/ { next } { getline seq } { print $0"\t"seq }' ${chimp_fasta_file} | sort -k1,1 |  awk '{print $0"\t"NR}'  > ${out_path}/${header_name}.chimp_coord.seq.table

human_fasta_file=${out_path}/${header_name}_human_coord_hg38.fa
awk '!/^>/ { next } { getline seq } { print $0"\t"seq }' ${human_fasta_file} | sort -k1,1 |  awk '{print $0"\t"NR}' > ${out_path}/${header_name}.human_coord.seq.table


paste -d "\t" <( awk '{print}' ${out_path}/${header_name}.chimp_coord.seq.table ) <( awk '{print}'  ${out_path}/${header_name}.human_coord.seq.table) > ${out_path}/${header_name}.seq.table.combined.temp

#check if the seq names match
awk '$1!=$4 {print}' ${out_path}/${header_name}.seq.table.combined.temp | head

awk '{print $1"\t"$2"\t"$5"\t"$3}' ${out_path}/${header_name}.seq.table.combined.temp > ${out_path}/${header_name}.seq.table.combined


#v2 is mostly changed to use the header name only and a seq table with both human and chimp sequence
#the first version is still usable
#submit and wait for run to finish
#bash ${script_dir}/get_mismatch_stats_v2.sh ${out_path} ${temp_dir} ${script_dir} ${header_name} ${out_path}/${header_name}.seq.table.combined
num_tasks_per_job=5000

#split_jobs divides up script so that each job runs ${num_tasks_per_job} tasks
bash ${script_dir}/get_mismatch_stats_v2_split_jobs_submit.sh ${out_path} ${temp_dir} ${script_dir} ${header_name} ${out_path}/${header_name}.seq.table.combined ${num_tasks_per_job}

################ merge files afterwards ###############

# merge all files afterwards 
out_file=${out_path}/${header_name}.human_coord.all_stats.txt

#bash ${script_dir}/merge_mismatch_stat_files_blast_v2.sh ${header_name} ${out_path} ${out_file} ${temp_dir} ${script_dir} ${out_path}/${header_name}.seq.table.combined

#simply cats files together
qsub -o ${temp_dir}/merge_mismatch_stats_submit_${header_name}.log ${script_dir}/merge_mismatch_stat_files_blast_v2.sh ${header_name} ${out_path} ${out_file} ${temp_dir} ${script_dir} ${out_path}/${header_name}.seq.table.combined
wc -l ${out_path}/${header_name}.human_coord.all_stats.txt

#add on seq length, as well as chromosome sampled from
inputFile=${out_path}/${header_name}.human_coord.all_stats.txt
paste <(awk '{print $1}' ${inputFile} ) <( awk '{print $1}' ${inputFile} | awk -F'[|]' '{print $4}' | awk -F'[_]' '{print $1}' ) <( awk -F"[\t|]" '{print $3-$2}' ${inputFile} ) <( cut -f2- ${inputFile} ) | sort -k1,1 > ${out_path}/${header_name}.human_coord.all_stats.revised.txt

################### get GC content proportions of MAF alignments - pulled from analyze_chip_LCL_data_v2.sh #######

####### this only needs to be run once #######################

mafdir=/cluster_path_temp/roadmap/LCL_data/all_maf_files
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
rm -fr ${mafdir}
mkdir -p ${mafdir}
refdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way
#${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt comes from finalFilteringProtocol_v4.sh
#use ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt because want positions where small human-specific deletions could fall

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}'  ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt > ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_chimp_cons_ids.txt
mafsInRegion -outDir ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_chimp_cons_ids.txt ${mafdir} ${refdir}/target_panTro4_multiz.maf

rm -f ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_chimp_cons_ids.txt
#output is number of sequences in each category, take the min.
#python ${del_dir}/getMAFStats.py ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt ${mafdir} ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats.txt
rm -f ${out_path}/getMAFStatsSubmit.log
qsub -o ${out_path}/getMAFStatsSubmit.log ${del_dir}/getMAFStatsSubmit.sh ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt ${mafdir} ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats.txt

#rearrange file to get conserved id only, with the statistics followed afterwards
#output columns are id, A%, T%, G%, C%, GC%, num in category, id showing all species in category, then repeated for all other categories, and followed by numMAFGAPS - which is actually zero for everything  
paste <(awk '{print $1"|"$2"|"$3"|"$4}' ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats.txt) <( cut -f5-  ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats.txt) > ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids.txt

#keep the relevant statistics of interest
#the initial statistics are the GC content for all blocks containing all 5 primates, and the number of species in the 5 primates
#awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"}' ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids.txt >  ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset.txt
#python checkDeletionOutputV3.py ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt ${mafdir} ${outputFile} ${windowSize}

#output columns are id, GC % for panTro4|rheMac8, num species in panTro4|rheMac8,  GC % for all 5 primates, num species in all 5 primates, %A, %T, %G, %C for all 5 primates, log odds score
awk '{print $1"\t"$6"\t"$7"\t"$20"\t"$21"\t"$16"\t"$17"\t"$18"\t"$19}' ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids.txt >  ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset.txt


##### add on info for log odds score from phastCons #########

#merge to get the log odds score?, first plot the log odds score against block length to see if its necessary
phastconsdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/cons/run.phast
awk '{print $1"|"$2"|"$3"\t"$5}' ${phastconsdir}/mostConserved.bed | sort -k1,1 > ${out_path}/mostConserved_additional_stats.txt
inputFile=${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset.txt
tempFile=${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset_temp.txt
outputFile=${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset_out.txt
paste <(awk '{print $1}'  ${inputFile} | awk -F"[|]" '{print $1"|"$2"|"$3}')  <(awk '{print}' ${inputFile} ) | sort -k1,1 > ${tempFile}
join -t $'\t' -11 ${tempFile} -21 ${out_path}/mostConserved_additional_stats.txt | cut -f2- > ${outputFile}
mv ${outputFile} ${inputFile}
rm -f ${tempFile}

########### get the same mismatch stats for the sequences with the deletions overlying them ####################
#I can possibly just do a join to get the stats?

#join to combine the conserved mapping to human statistics with the seq length statistics
#use the following file for delGapStats.R?, an initial metadata with files of interest
join -t$'\t' -11 -21 <(sort -k1,1 ${out_path}/${header_name}.human_coord.all_stats.revised.txt) <( sort -k1,1 ${out_path}/mostConservedMapped_new_0.001_human_size_filtered_cons_stats_with_ids_subset.txt ) > ${out_path}/${header_name}_all_stats.txt


############### join to get deletions tested/cons stats ##################

awk 'NR%2==1{print}' ${fasta_file} | grep -E -v "EMVAR|LCL_HEPG2_POS" | awk -F"[>|]" '{print $8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14"|"$15"|"$16"|"$17"|"$18"|"$19"|"$20"|"$21}' | sort | uniq > ${out_path}/${real_data_header}_id.txt
awk -F"|" '{print $1"|"$2"|"$3"|"$4"|"$5"|"$6}' ${out_path}/${real_data_header}_id.txt | sort | uniq > ${out_path}/${real_data_header}.orig_cons_ID.txt

join -t$'\t' -11 -21 ${out_path}/${real_data_header}.orig_cons_ID.txt ${out_path}/${header_name}_all_stats.txt  > ${out_path}/${real_data_header}.cons_stats.txt

#*********** rename files if needed ***********#

#remove individual files

qsub -o ${temp_dir}/remove_intermediate_mismatch_blast_files.log ${script_dir}/remove_intermediate_mismatch_blast_files.sh ${perm_header_name} ${out_path} ${out_path}/${perm_header_name}.seq.table.combined

#change old names if needed

ls -l ${out_path}/*${perm_header_name}* | awk '{print $NF}' > ${out_path}/replace_names_files_1.txt

replace_name=${mostConservedMapped_new_0.001_human_size_filtered_extended_${maxExtendConsSeqLen}_hg38_aligned_panTro4_check_stats_revised}
sed "s#${perm_header_name}#${header_name}#g" ${out_path}/replace_names_files_1.txt > ${out_path}/replace_names_files_2.txt

paste -d "\t" ${out_path}/replace_names_files_1.txt ${out_path}/replace_names_files_2.txt > ${out_path}/replace_names_files.txt
while read -a array
do
	orig_file_name="${array[0]}"
	new_file_name="${array[1]}"
	mv ${orig_file_name} ${new_file_name}
done<${out_path}/replace_names_files.txt
