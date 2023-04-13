#****************** preprocess permutation files for intersection - no need to run if already ran ******************#
#this parameter is for script extendCoordinatesOligoDesignConsBackground.py 
maxExtendConsSeqLen=200
#extend coordinates all coordinates to 200 bp length, if over ${maxExtendConsSeqLen} don't extend
#we use this parameter to align the human sequence to the chimp sequence as a control for number of mismatches

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

perm_type_results_dir=/cluster_path/ape_project/deletions_project/permutation_type_results

#parse perm human hg19 coord, create file containing all the permutation files
permutations_file_path=/cluster_path/ape_project/deletions_project/roadmap/permutations_new
perm_seq_header="mostCons_all_mismatch_analysis_maxExtendLen_${maxExtendConsSeqLen}_lengthCf_${lengthBoundPct}_mismatchCf_${misMatchBoundPct}_GCCf_${GCBoundPct}_logOddsCf_${logOddsBoundPct}_maxTrys_${sampleNullSeqMaxTrys}_devPct_${sampleNullSeqConsDevPct}"
numPerm=1000

rm -f ${perm_type_results_dir}/hcondel_general_category_human_hg19_all_input_perm_files.txt
for permIter in $(seq 1 ${numPerm})
do
    echo -e "${permutations_file_path}/${permIter}_${perm_seq_header}_list.txt" >> ${perm_type_results_dir}/hcondel_general_category_human_hg19_all_input_perm_files.txt
done

#parse perm chimp panTro4 coord, create file containing all the permutation files
permutations_file_path=/cluster_path/ape_project/deletions_project/tf_enrichment/files/permutations
perm_seq_header="mostCons_all_mismatch_analysis_maxExtendLen_${maxExtendConsSeqLen}_lengthCf_${lengthBoundPct}_mismatchCf_${misMatchBoundPct}_GCCf_${GCBoundPct}_logOddsCf_${logOddsBoundPct}_maxTrys_${sampleNullSeqMaxTrys}_devPct_${sampleNullSeqConsDevPct}"
numPerm=1000

rm -f ${perm_type_results_dir}/hcondel_general_category_chimp_panTro4_all_input_perm_files.txt
for permIter in $(seq 1 ${numPerm})
do
    #rm -f ${permutations_file_path}/${permIter}_${perm_seq_header}_chimp_coord.no_strand.txt
    #awk '{print $1"\t"$2"\t"$3"\t"$5}' ${permutations_file_path}/${permIter}_${perm_seq_header}_chimp_coord.txt | sort -k1,1 -k2,2n > ${permutations_file_path}/${permIter}_${perm_seq_header}_chimp_coord.no_strand.txt
    echo -e "${permutations_file_path}/${permIter}_${perm_seq_header}_chimp_coord.no_strand.txt" >> ${perm_type_results_dir}/hcondel_general_category_chimp_panTro4_all_input_perm_files.txt
done

#parse perm human hg38 coord, create file containing all the permutation files
permutations_file_path=/cluster_path/ape_project/deletions_project/tf_enrichment/files/permutations
perm_seq_header="mostCons_all_mismatch_analysis_maxExtendLen_${maxExtendConsSeqLen}_lengthCf_${lengthBoundPct}_mismatchCf_${misMatchBoundPct}_GCCf_${GCBoundPct}_logOddsCf_${logOddsBoundPct}_maxTrys_${sampleNullSeqMaxTrys}_devPct_${sampleNullSeqConsDevPct}"
numPerm=1000

rm -f ${perm_type_results_dir}/hcondel_general_category_human_hg38_all_input_perm_files.txt
for permIter in $(seq 1 ${numPerm})
do
    rm -f ${permutations_file_path}/${permIter}_${perm_seq_header}_human_coord.no_strand.txt
    awk '{print $1"\t"$2"\t"$3"\t"$5}' ${permutations_file_path}/${permIter}_${perm_seq_header}_human_coord.txt | sort -k1,1 -k2,2n > ${permutations_file_path}/${permIter}_${perm_seq_header}_human_coord.no_strand.txt
    echo -e "${permutations_file_path}/${permIter}_${perm_seq_header}_human_coord.no_strand.txt" >> ${perm_type_results_dir}/hcondel_general_category_human_hg38_all_input_perm_files.txt
done


#********************************* rest of analysis ************************************************#
script_dir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

#for age, if multiple blocks/get oldest age

del_dir=/cluster_path/ape_project/deletions_project
age_dir=${del_dir}/Ages

repeat_dir=/cluster_path/ape_project/deletions_project/repeats

screen_ccre_dir=/cluster_path/ape_project/deletions_project/screen_ccre

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1

perm_type_results_dir=/cluster_path/ape_project/deletions_project/permutation_type_results
mkdir -p ${perm_type_results_dir}

#need to manually make ${perm_type_results_dir}/type_info.age_repeat.txt
#type_ref_file=${perm_type_results_dir}/type_info.age_repeat.txt

#cat ${perm_type_results_dir}/type_info.age_repeat.txt > ${perm_type_results_dir}/type_info.txt

#cat ${perm_type_results_dir}/type_info.genomic_regions.txt >> ${perm_type_results_dir}/type_info.txt

#gtex v8
#cat ${perm_type_results_dir}/type_info.gtex_v8.txt >> ${perm_type_results_dir}/type_info.txt

#for all
type_ref_file=${perm_type_results_dir}/type_info.txt

##### construct meta file for just screen ccres #####
awk '{print $2}' ${screen_ccre_dir}/all_ccre_info.txt > ${screen_ccre_dir}/all_ccre_info.screen_ccre_file_labels.txt

#hg19 sorted bed file
awk '{print $1}' ${screen_ccre_dir}/all_ccre_info.txt | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 1; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}' > ${screen_ccre_dir}/all_ccre_info.orig_header_names.txt
awk -v screen_ccre_dir=${screen_ccre_dir} '{print screen_ccre_dir"/"$1".parsed.hg19_mapped.sorted.bed"}' ${screen_ccre_dir}/all_ccre_info.orig_header_names.txt > ${screen_ccre_dir}/all_ccre_info.hg19_sorted_bed_files.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$2"_labels_uniq.txt"}' ${screen_ccre_dir}/all_ccre_info.txt >  ${screen_ccre_dir}/all_ccre_info.labels_file_list.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$1"_perm_initial_results_file_list.txt"}' ${screen_ccre_dir}/all_ccre_info.screen_ccre_file_labels.txt >  ${screen_ccre_dir}/all_ccre_info.perm_file_list.txt

num_lines=$(wc -l ${screen_ccre_dir}/all_ccre_info.screen_ccre_file_labels.txt | awk '{print $1}')
paste <(cat ${screen_ccre_dir}/all_ccre_info.screen_ccre_file_labels.txt) <(cat ${screen_ccre_dir}/all_ccre_info.hg19_sorted_bed_files.txt) <(cat ${screen_ccre_dir}/all_ccre_info.labels_file_list.txt) <(cat ${screen_ccre_dir}/all_ccre_info.perm_file_list.txt) <(for (( c=1; c<=num_lines; c++)) ; do echo -e "${perm_type_results_dir}/hcondel_general_category_human_hg19_all_input_perm_files.txt"; done) > ${perm_type_results_dir}/type_info.screen_ccre.txt 

##### construct meta file for just genomic regions #####
genomic_regions_dir=/cluster_path/ape_project/deletions_project/genomic_regions

awk '{print $2}' ${genomic_regions_dir}/all_genomic_regions_info.txt > ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions_file_labels.txt

#sorted bed file
awk -v genomic_regions_dir=${genomic_regions_dir} '{print genomic_regions_dir"/"$1".sorted.bed"}' ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions_file_labels.txt  > ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions.sorted_bed_files.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$2"_labels_uniq.txt"}' ${genomic_regions_dir}/all_genomic_regions_info.txt >  ${genomic_regions_dir}/all_genomic_regions_info.labels_file_list.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$1"_perm_initial_results_file_list.txt"}' ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions_file_labels.txt >  ${genomic_regions_dir}/all_genomic_regions_info.perm_file_list.txt

num_lines=$(wc -l ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions_file_labels.txt | awk '{print $1}')
paste <(cat ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions_file_labels.txt) <(cat ${genomic_regions_dir}/all_genomic_regions_info.genomic_regions.sorted_bed_files.txt) <(cat ${genomic_regions_dir}/all_genomic_regions_info.labels_file_list.txt) <(cat ${genomic_regions_dir}/all_genomic_regions_info.perm_file_list.txt) <(for (( c=1; c<=num_lines; c++)) ; do echo -e "${perm_type_results_dir}/hcondel_general_category_chimp_panTro4_all_input_perm_files.txt"; done) > ${perm_type_results_dir}/type_info.genomic_regions.txt 

##### construct meta file for just GTEx v8 #####
#see GTEx preprocess in code below for making ${gtex_v8_dir}/all_gtex_v8_info.txt
gtex_v8_dir=/cluster_path/ape_project/deletions_project/gtex_v8

awk '{print $2}' ${gtex_v8_dir}/all_gtex_v8_info.txt > ${gtex_v8_dir}/all_gtex_v8_info.file_labels.txt

#sorted bed file
awk -v gtex_v8_dir=${gtex_v8_dir} '{print gtex_v8_dir"/"$1".sorted.bed"}' ${gtex_v8_dir}/all_gtex_v8_info.file_labels.txt > ${gtex_v8_dir}/all_gtex_v8_info.sorted_bed_files.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$2"_labels_uniq.txt"}' ${gtex_v8_dir}/all_gtex_v8_info.txt > ${gtex_v8_dir}/all_gtex_v8_info.labels_file_list.txt

#perm file list
awk -v perm_type_results_dir=${perm_type_results_dir} '{print perm_type_results_dir"/"$1"_perm_initial_results_file_list.txt"}' ${gtex_v8_dir}/all_gtex_v8_info.file_labels.txt > ${gtex_v8_dir}/all_gtex_v8_info.perm_file_list.txt

num_lines=$(wc -l ${gtex_v8_dir}/all_gtex_v8_info.file_labels.txt | awk '{print $1}')
paste <(cat ${gtex_v8_dir}/all_gtex_v8_info.file_labels.txt) <(cat ${gtex_v8_dir}/all_gtex_v8_info.sorted_bed_files.txt) <(cat ${gtex_v8_dir}/all_gtex_v8_info.labels_file_list.txt) <(cat ${gtex_v8_dir}/all_gtex_v8_info.perm_file_list.txt) <(for (( c=1; c<=num_lines; c++)) ; do echo -e "${perm_type_results_dir}/hcondel_general_category_human_hg38_all_input_perm_files.txt"; done) > ${perm_type_results_dir}/type_info.gtex_v8.txt 

#make file for measuring different attributes of same label, this can be for example distance to overlap label
#modify this for different iterations
while read -a array
do
    type_iter="${array[0]}"
    label_ids_file="${array[2]}"

    num_lines=$(wc -l ${gtex_v8_dir}/gtex_v8_dist_cutoff_interest_file.txt | awk '{print $1}')
    while read -a array2
    do
        label="${array2[0]}"
        label_file_id="${array2[1]}"
        paste <(for (( c=1; c<=num_lines; c++)) ; do echo -e "${label}\t${c}\t${label_file_id}"; done) <(cat ${gtex_v8_dir}/gtex_v8_dist_cutoff_interest_file.txt) | awk '{print $1"\t"$2"\t"$3"_distance_"$4}' > "${label_ids_file}.count_additional_meta.${label_file_id}"
    done<${label_ids_file}
done<${perm_type_results_dir}/type_info.gtex_v8.txt 

#*************** preprocessing (can skip if already done) ****************************#

####### lift denylist filtered to hg19 ########

denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt
paste <(awk -F"|" '{print $10"\t"$11-1"\t"$12}' ${denylist_filtered_ids_file}) <(awk '{print $0}'  ${denylist_filtered_ids_file}) > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.txt

chain_file=/cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz
liftOver ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.txt ${chain_file} ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.txt ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_unmapped.txt

sort -k1,1 -k2,2n ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.txt > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.txt 
sort -k1,1 -k2,2n ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.txt > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.sorted.txt

####### panTro4 chimp coords ########

denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt
paste <(awk -F"|" '{print $7"\t"$8"\t"$9}' ${denylist_filtered_ids_file}) <(awk '{print $0}'  ${denylist_filtered_ids_file}) | sort -k1,1 -k2,2n > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_chimp_coord.sorted.txt

####### age_del preprocess ###########
#sort -k1,1 -k2,2n ${age_dir}/hg19_aged.txt > ${age_dir}/hg19_aged.sorted.txt
type_iter=age_del

#get all age labels
awk '{print $1}' ${age_dir}/sorted_ages.rev.txt| sort -k1,1 | uniq > ${perm_type_results_dir}/age_del_labels_uniq.txt
paste <(cat ${perm_type_results_dir}/age_del_labels_uniq.txt) <(cat ${perm_type_results_dir}/age_del_labels_uniq.txt) > ${perm_type_results_dir}/age_del_labels_uniq.temp.txt
mv ${perm_type_results_dir}/age_del_labels_uniq.temp.txt ${perm_type_results_dir}/age_del_labels_uniq.txt

type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
join -11 <(sort -k1,1 ${denylist_filtered_ids_file}) -21 <(sort -k1,1 ${age_dir}/all_files_af_${af_filter}_final_filtered.del_age_info.txt) | awk '{print $1"\t"$2"\t"$4}' | sort -k1,1 -k2,2n > ${perm_type_results_dir}/${type_header_id}_info.txt

####### repeat_class preprocess ###########

type_iter=repeat_class
type_bed_file=${repeat_dir}/rpt_withClass.hg19_mapped.sorted.txt

#map to hg19
chain_file=/cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz
liftOver ${repeat_dir}/rpt_withClass.sorted.txt ${chain_file} ${repeat_dir}/rpt_withClass.hg19_mapped.txt ${repeat_dir}/rpt_withClass.hg19_unmapped.txt

sort -k1,1 -k2,2n ${repeat_dir}/rpt_withClass.hg19_mapped.txt > ${repeat_dir}/rpt_withClass.hg19_mapped.sorted.txt

awk '{print $NF}' ${repeat_dir}/rpt_withClass.hg19_mapped.sorted.txt | sort -k1,1 | uniq > ${perm_type_results_dir}/repeat_class_labels_uniq.txt
paste <(cat ${perm_type_results_dir}/repeat_class_labels_uniq.txt) <(cat ${perm_type_results_dir}/repeat_class_labels_uniq.txt | sed 's#/#_#g') > ${perm_type_results_dir}/repeat_class_labels_uniq.temp.txt
mv ${perm_type_results_dir}/repeat_class_labels_uniq.temp.txt ${perm_type_results_dir}/repeat_class_labels_uniq.txt

#intersect with hcondel file
type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
posFile=${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.txt 
bedtools intersect -sorted -a ${posFile} -b ${type_bed_file} -loj > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.${type_iter}_intersected.txt

python ${del_dir}/get_random_intersected_label.py ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.${type_iter}_intersected.txt ${perm_type_results_dir}/${type_header_id}_info.txt

####### screen_ccre preprocess ###########
all_regions_info_file=/cluster_path/ape_project/deletions_project/screen_ccre/all_ccre_info.txt
#for all screen ccres
while read -a array
do
	download_link=${array[0]}
	type_iter=${array[1]}
	
	orig_file_header=$(echo ${download_link} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 1; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}')
	
	type_bed_file=${screen_ccre_dir}/${orig_file_header}.parsed.hg19_mapped.sorted.bed
	
	awk '{print $1"\t"$2"\t"$3"\t1"}' ${screen_ccre_dir}/${orig_file_header}.bed  > ${screen_ccre_dir}/${orig_file_header}.parsed.bed 
	chain_file=/cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz
	liftOver ${screen_ccre_dir}/${orig_file_header}.parsed.bed ${chain_file} ${screen_ccre_dir}/${orig_file_header}.parsed.hg19_mapped.bed ${screen_ccre_dir}/${orig_file_header}.parsed.hg19_unmapped.bed 

	sort -k1,1 -k2,2n ${screen_ccre_dir}/${orig_file_header}.parsed.hg19_mapped.bed > ${screen_ccre_dir}/${orig_file_header}.parsed.hg19_mapped.sorted.bed
	echo -e "1\t${type_iter}" > ${perm_type_results_dir}/${type_iter}_labels_uniq.txt

	#intersect with hcondel file
	type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
	posFile=${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.txt 
	bedtools intersect -sorted -a ${posFile} -b ${type_bed_file} -loj > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.${type_iter}_intersected.txt

	python ${del_dir}/get_random_intersected_label.py ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.hg19_mapped.sorted.${type_iter}_intersected.txt ${perm_type_results_dir}/${type_header_id}_info.txt

	#cp ${script_dir}/get_perm_stats_screen_ccre_foundation.sh ${script_dir}/get_perm_stats_${orig_file_header}.sh 
done<${all_regions_info_file}


#for all screen ccres, copy script, and modify to adapt
while read -a array
do
	download_link=${array[0]}
	type_iter=${array[1]}
	
	orig_file_header=$(echo ${download_link} | awk -F"/" '{print $NF}' | awk -F"." '{for (i = 1; i <= NF-2; i++) {printf $i"."}; printf $(NF-1)"\n"}')
	
	cp ${script_dir}/get_perm_stats_screen_ccre_foundation.sh ${script_dir}/get_perm_stats_${type_iter}.sh 
	
	sed -i "s#screen_ccre#$type_iter#g" ${script_dir}/get_perm_stats_${type_iter}.sh 
done<${all_regions_info_file}

####### genomic_regions preprocess ###########
all_regions_info_file=/cluster_path/ape_project/deletions_project/genomic_regions/all_genomic_regions_info.txt
#for all genomic regions
while read -a array
do
	initial_bed_file=${array[0]}
	type_iter=${array[1]}
	
	sort -k1,1 -k2,2n ${initial_bed_file} | awk -v type_iter=${type_iter} '{print $1"\t"$2"\t"$3"\t1"}' > ${genomic_regions_dir}/${type_iter}.sorted.bed
	
	type_bed_file=${genomic_regions_dir}/${type_iter}.sorted.bed
	echo -e "1\t${type_iter}" > ${perm_type_results_dir}/${type_iter}_labels_uniq.txt

	#intersect with hcondel file
	type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
	posFile=${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_chimp_coord.sorted.txt
	bedtools intersect -sorted -a ${posFile} -b ${type_bed_file} -loj > ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_chimp_coord.sorted.${type_iter}_intersected.txt

	python ${del_dir}/get_random_intersected_label.py ${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_chimp_coord.sorted.${type_iter}_intersected.txt ${perm_type_results_dir}/${type_header_id}_info.txt

	#cp ${script_dir}/get_perm_stats_genomic_regions_foundation.sh ${script_dir}/get_perm_stats_${orig_file_header}.sh 
done<${all_regions_info_file}


#for all genomic regions, copy script, and modify to adapt
while read -a array
do
	type_iter=${array[1]}
	
	cp ${script_dir}/get_perm_stats_screen_ccre_foundation.sh ${script_dir}/get_perm_stats_${type_iter}.sh 
	
	sed -i "s#screen_ccre#$type_iter#g" ${script_dir}/get_perm_stats_${type_iter}.sh 
done<${all_regions_info_file}

####### GTEx preprocess ###########

### initial parameters ###
gtex_v8_dir=/cluster_path/ape_project/deletions_project/gtex_v8
del_dir=/cluster_path/ape_project/deletions_project
gtex_v8_temp_dir=/cluster_path_temp/gtex_v8
mkdir -p ${gtex_v8_dir}
mkdir -p ${gtex_v8_dir}/scripts
mkdir -p ${gtex_v8_dir}/output
mkdir -p ${gtex_v8_temp_dir}

#download two files below from GTEx, see below for download
sample_gene_exp_table_file=${gtex_v8_temp_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
sample_attributes_table_file=${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.count_table_intersect.txt

all_broad_tissue_file=${gtex_v8_dir}/GTEx_Analysis_v8_all_broad_tissues.txt
all_specific_tissue_file=${gtex_v8_dir}/GTEx_Analysis_v8_all_specific_tissues.txt

fdr_log2FC_cutoff_info_file_broad=${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_broad.txt
fdr_log2FC_cutoff_info_file_specific=${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_specific.txt

rm -f ${fdr_log2FC_cutoff_info_file_broad}
rm -f ${fdr_log2FC_cutoff_info_file_specific}

echo -e "0.1\t10" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t9" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t8" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t7" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t6" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t5" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t4" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t3" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t2" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t1" >> ${fdr_log2FC_cutoff_info_file_broad}
echo -e "0.1\t0.5" >> ${fdr_log2FC_cutoff_info_file_broad}

cp ${fdr_log2FC_cutoff_info_file_broad} ${fdr_log2FC_cutoff_info_file_specific}

dist_cutoff_interest_file=${gtex_v8_dir}/gtex_v8_dist_cutoff_interest_file.txt
rm -f ${dist_cutoff_interest_file}

echo -e "50" >> ${dist_cutoff_interest_file}
echo -e "100" >> ${dist_cutoff_interest_file}
echo -e "500" >> ${dist_cutoff_interest_file}
echo -e "1000" >> ${dist_cutoff_interest_file}
echo -e "5000" >> ${dist_cutoff_interest_file}
echo -e "10000" >> ${dist_cutoff_interest_file}
echo -e "50000" >> ${dist_cutoff_interest_file}
echo -e "100000" >> ${dist_cutoff_interest_file}
echo -e "500000" >> ${dist_cutoff_interest_file}


#download gene tpm matrix
#wget -P ${gtex_v8_dir} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
#gunzip -c ${gtex_v8_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz > ${gtex_v8_temp_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct

#download gene count matrix
wget -P ${gtex_v8_dir} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
gunzip -c ${gtex_v8_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz > ${gtex_v8_temp_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct

#download sample annotation matrix
wget -P ${gtex_v8_dir} https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

#download gtf file 
wget -P ${gtex_v8_dir} https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf

#get coordinates from gtf file
awk '{if($3=="gene" && $0 ~/gene_type \"protein_coding\"/){print $1"\t"$4-1"\t"$5"\t"$7"\t"$10"\t"$16}}' ${gtex_v8_dir}/gencode.v26.GRCh38.genes.gtf | tr -d '";' > ${gtex_v8_dir}/gencode.v26.GRCh38.genes.gtf.protein_coding.bed

#keep only tissues with columns in the gene count matrix
head -3 ${gtex_v8_temp_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct | awk 'NR==3{print}' | tr -s ' \t' '\n' | awk 'NR>2{print}' | sort -t$'\t' -k1,1 > ${gtex_v8_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.colnames
awk -F"\t" 'NR>1{print $1"\t"$6"\t"$7}' ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | sort -t$'\t' -k1,1 > ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.txt
join -t $'\t' -11 <(cat ${gtex_v8_dir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.colnames) -21 <(cat ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.txt) > ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.count_table_intersect.txt

#get file names of all broad tissues
awk -F"\t" 'NR>1{print $2}' ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.count_table_intersect.txt | sort -t$'\t' | uniq > ${gtex_v8_dir}/GTEx_Analysis_v8_all_broad_tissues.txt

#get file names of all specific tissues, keep only specific tissues that has at least two other specific tissues in larger domain
awk -F"\t" 'NR>1{print $2"\t"$3}' ${gtex_v8_dir}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.count_table_intersect.txt | sort -t$'\t' | uniq > ${gtex_v8_dir}/GTEx_v8_tissue_map.txt
awk -F"\t" '{print $1}' ${gtex_v8_dir}/GTEx_v8_tissue_map.txt | sort -t$'\t' | uniq -c | sed -E 's/^ *//; s/ /\t/'| awk -F"\t" '$1>1{print $2}' > ${gtex_v8_dir}/GTEx_v8_tissue_map.mult_specific_tissue_list.txt
join -t $'\t' -11 -21 <(cat ${gtex_v8_dir}/GTEx_v8_tissue_map.mult_specific_tissue_list.txt) <(sort -t$'\t' -k1,1 ${gtex_v8_dir}/GTEx_v8_tissue_map.txt) | awk -F"\t" '{print $2}' >  ${gtex_v8_dir}/GTEx_Analysis_v8_all_specific_tissues.txt
rm -f ${gtex_v8_dir}/GTEx_v8_tissue_map.mult_specific_tissue_list.txt

##

#get pvals/log2FC of all broad tissues
out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_broad_tissues
bash ${gtex_v8_dir}/scripts/get_preferred_exp_all_tissues_submit.sh ${sample_gene_exp_table_file} ${sample_attributes_table_file} ${out_file_path_with_file_header} ${all_broad_tissue_file}

#get pvals/log2FC of specific tissues from larger domains
out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_specific_tissues
bash ${gtex_v8_dir}/scripts/get_preferred_exp_all_tissues_submit.sh ${sample_gene_exp_table_file} ${sample_attributes_table_file} ${out_file_path_with_file_header} ${all_specific_tissue_file}

#head ${out_file_path_with_file_header}_get_preferred_exp.log
#head ${out_file_path_with_file_header}_get_preferred_exp_all_runs.txt

##

#filter out GTEx based off of cutoffs

#get TSS
awk '{if($4=="+"){print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} else{print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}}' ${gtex_v8_dir}/gencode.v26.GRCh38.genes.gtf.protein_coding.bed > ${gtex_v8_dir}/gencode.v26.GRCh38.genes.gtf.protein_coding.TSS.bed
gene_meta_table_file=${gtex_v8_dir}/gencode.v26.GRCh38.genes.gtf.protein_coding.TSS.bed

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_broad_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array
    do
        tissue_of_interest="${array[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        join -t $'\t' -15 -21 <(sort -k5,5 ${gene_meta_table_file}) <(awk -F"\t" 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' ${out_file_path_with_file_header}_${tissue_of_interest}.txt | sort -k1,1) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9}' > ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt

        #filter out genes based on cutoffs
        awk -F"\t" -v fdr_cutoff="${fdr_cutoff}" -v log2FC_cutoff="${log2FC_cutoff}" '($7!="NA") && ($9!="NA") && ($9<fdr_cutoff) && ($7>log2FC_cutoff || $7< -log2FC_cutoff) {print}' ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt > ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}.txt

    done<${all_broad_tissue_file}
done<${fdr_log2FC_cutoff_info_file_broad}

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_specific_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array
    do
        tissue_of_interest="${array[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        join -t $'\t' -15 -21 <(sort -k5,5 ${gene_meta_table_file}) <(awk -F"\t" 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' ${out_file_path_with_file_header}_${tissue_of_interest}.txt | sort -k1,1) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9}' > ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt

        #filter out genes based on cutoffs
        awk -F"\t" -v fdr_cutoff="${fdr_cutoff}" -v log2FC_cutoff="${log2FC_cutoff}" '($7!="NA") && ($9!="NA") && ($9<fdr_cutoff) && ($7>log2FC_cutoff || $7< -log2FC_cutoff) {print}' ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt > ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}.txt

    done<${all_specific_tissue_file}
done<${fdr_log2FC_cutoff_info_file_specific}


## create file for permutation overlap ##
all_regions_info_file=/cluster_path/ape_project/deletions_project/gtex_v8/all_gtex_v8_info.txt
rm -f ${all_regions_info_file}

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_broad_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array
    do
        tissue_of_interest="${array[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        out_iter_file="${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}.txt"
    
        #check if any genes exists
        num_sig_genes=$(wc -l ${out_iter_file} | awk '{print $1}')
        if [ ${num_sig_genes} -gt 0 ]
        then
            echo -e "${out_iter_file}\tgtex_v8_${tissue_of_interest}_fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}" >> ${all_regions_info_file}
        else
            echo -e "No significant genes in ${tissue_of_interest}"
        fi
    
    done<${all_broad_tissue_file}
done<${fdr_log2FC_cutoff_info_file_broad}

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_specific_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array
    do
        tissue_of_interest="${array[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        out_iter_file="${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}.txt"
    
        #check if any genes exists
        num_sig_genes=$(wc -l ${out_iter_file} | awk '{print $1}')
        if [ ${num_sig_genes} -gt 0 ]
        then
            echo -e "${out_iter_file}\tgtex_v8_${tissue_of_interest}_fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}" >> ${all_regions_info_file}
        else
            echo -e "No significant genes in ${tissue_of_interest}"
        fi
    
    done<${all_specific_tissue_file}
done<${fdr_log2FC_cutoff_info_file_specific}

all_regions_info_file=/cluster_path/ape_project/deletions_project/gtex_v8/all_gtex_v8_info.txt
#for all tissues
while IFS=$'\t' read -a array
do
    initial_bed_file="${array[0]}"
    type_iter="${array[1]}"

    sort -k1,1 -k2,2n ${initial_bed_file} | awk -v type_iter=${type_iter} '{print $1"\t"$2"\t"$3"\t1"}' > ${gtex_v8_dir}/${type_iter}.sorted.bed

    type_bed_file=${gtex_v8_dir}/${type_iter}.sorted.bed
    echo -e "1\t${type_iter}" > ${perm_type_results_dir}/${type_iter}_labels_uniq.txt

    #intersect with hcondel file
    type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
    posFile=${perm_type_results_dir}/novaseq_12_15_18_denylist_filtered_human_coord.sorted.txt

    #get closest feature, parse to see if it passes distance cutoffs in ${dist_cutoff_interest_file} 
    bedtools closest -d -a ${posFile} -b ${type_bed_file} | awk '{print $4"\t"$NF}' > "${perm_type_results_dir}/${type_header_id}_info.temp.txt"
    python ${script_dir}/format_bedtools_closest_dist_output.py "${perm_type_results_dir}/${type_header_id}_info.temp.txt" "${dist_cutoff_interest_file}" "${perm_type_results_dir}/${type_header_id}_info.txt"
    rm -f "${perm_type_results_dir}/${type_header_id}_info.temp.txt"
done<${all_regions_info_file}


#for all tissues, copy script, and modify to adapt
while IFS=$'\t' read -a array
do
    type_iter="${array[1]}"

    cp ${script_dir}/get_perm_stats_closest_dist_foundation_v2.sh ${script_dir}/get_perm_stats_${type_iter}.sh 

    sed -i "s#sample#$type_iter#g; s#type#gtex_v8#g" ${script_dir}/get_perm_stats_${type_iter}.sh 
done<${all_regions_info_file}

##################### get perm stats ################

#wrapper script to submit and get stats from all permutation runs
rm -f ${script_dir}/submit_all_get_perm_stats_type_runs.log 
qsub -o ${script_dir}/submit_all_get_perm_stats_type_runs.log ${script_dir}/submit_all_get_perm_stats_type_runs.sh ${type_ref_file} ${perm_type_results_dir} ${script_dir}

#after running above, combine files
bash ${script_dir}/combine_perm_stats_type_runs_submit.sh ${type_ref_file} ${perm_type_results_dir} ${script_dir}
#head ${script_dir}/combine_perm_stats_type_runs.log

### get count stats of permutations across all classes ###

## script to run above

#nohup bash ${script_dir}/get_perm_count_stats.sh &> ${perm_type_results_dir}/score_summary.out&

grep -v "gtex_v8" ${perm_type_results_dir}/type_info.txt > ${type_ref_file}.count
type_ref_file_count=${type_ref_file}.count

bash ${script_dir}/get_perm_count_stats_submit.sh ${type_ref_file_count} ${perm_type_results_dir} ${script_dir}
#head ${script_dir}/get_perm_count_stats_all.log
#${script_dir}/get_perm_count_stats.sh

### compare with empirical, output summary file ###

type_ref_file=${perm_type_results_dir}/type_info.txt
bash ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_submit.sh ${type_ref_file_count} ${perm_type_results_dir} ${script_dir}
tail ${script_dir}/perm_analyze_stats_general_count_enrichment_v2_all.log

## make final file  ##
type_ref_file=${perm_type_results_dir}/type_info.txt
rm -f ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt
echo -e "label\tempirical_count\tcount_pval\tperm_count_mean\tperm_count_sd\tcount_z_score\tclass" > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt
while read -a array
do
    type_iter=${array[0]}
    out_path=${perm_type_results_dir}/${type_iter}
    cat ${out_path}/${type_iter}_count_enrichment_hcondel_denylist_filtered_final_stats.txt >> ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt
done<${type_ref_file_count}

### get score stats ###
#the last field contains any score

grep -E "gtex_v8" ${type_ref_file} > ${type_ref_file}.score
type_ref_file_score=${type_ref_file}.score

bash ${script_dir}/get_perm_score_stats_submit.sh ${type_ref_file_score} ${perm_type_results_dir} ${script_dir}

#after running, combine output score files

outFile=${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.txt

rm -f ${outFile}
echo -e "type\tempirical_score_avg\tperm_score_avg_pval\tperm_score_avg_mean\tperm_score_avg_sd\tperm_score_avg_z_score" > ${outFile}
while read -a array
do
    type_iter=${array[0]}

    type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_empirical
    real_stats_file=${perm_type_results_dir}/${type_header_id}_info.txt

    type_header_id=${type_iter}_enrichment_hcondel_denylist_filtered_random
    perm_summary_stats_file=${perm_type_results_dir}/${type_iter}/${type_header_id}_score_summary_stats.txt

    if [[ ${type_iter} =~ "gtex_v8" ]]; then
        real_score_avg=$(awk '$NF>-1{ sum += $NF; n++ } END { if (n > 0) print sum / n; }' ${real_stats_file})
    else
        real_score_avg=$(awk '{ sum += $NF; n++ } END { if (n > 0) print sum / n; }' ${real_stats_file})
    fi

    perm_score_avg_num=$(awk -F"\t" -v real_score_avg=${real_score_avg} '$1>real_score_avg{print}' ${perm_summary_stats_file} | wc -l)
    perm_score_avg_denom=$(awk '{print}' ${perm_summary_stats_file} | wc -l)

    if [ ${perm_score_avg_denom} -eq 0 ]
    then
        perm_score_avg_pval="NA"
    else
        perm_score_avg_pval=$(echo -e "${perm_score_avg_num}\t${perm_score_avg_denom}" | awk -F"\t" '{print $1/$2}')
    fi
    
    perm_score_avg_mean=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${perm_summary_stats_file})
    perm_score_avg_sd=$(awk '{sum += $1; sumsq += ($1)^2} END {print sqrt((sumsq-sum^2/NR)/(NR-1))}' ${perm_summary_stats_file})

    #https://stackoverflow.com/questions/2424770/floating-point-comparison-in-shell
    var=$(awk 'BEGIN{ print "'${perm_score_avg_sd}'"==0 }')
    if [ ${var} -eq 1 ]
    then
        perm_score_avg_z_score="NA"
    else
        perm_score_avg_z_score=$( echo -e "${real_score_avg}\t${perm_score_avg_mean})\t${perm_score_avg_sd}" | awk -F"\t" '{print ($1-$2)/$3}' )
    fi

    
    echo -e "${type_iter}\t${real_score_avg}\t${perm_score_avg_pval}\t${perm_score_avg_mean}\t${perm_score_avg_sd}\t${perm_score_avg_z_score}" >> ${outFile}

done<${type_ref_file_score}	

#rm -f ${script_dir}/get_score_summary.log 
#qsub -o ${script_dir}/get_score_summary.log ${script_dir}/get_score_summary.sh

#look at data

awk 'NR>1{print}' ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.txt | grep gtex | sort -k3,3n | awk '$3>0.95 {print}'
awk 'NR>1{print}' ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt | grep gtex | sort -k3,3n | awk '$3<0.05 {print}'

awk 'NR>1{print}' ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt | grep gtex | sort -k3,3n | awk '$3<0.05 && $6!="NA" && $4>5 {print}'

####### for gtex_v8, create master summary file including the number of genes in each fdr and log2FC cutoff #######

gtex_v8_dir=/cluster_path/ape_project/deletions_project/gtex_v8
all_broad_tissue_file=${gtex_v8_dir}/GTEx_Analysis_v8_all_broad_tissues.txt
all_specific_tissue_file=${gtex_v8_dir}/GTEx_Analysis_v8_all_specific_tissues.txt

fdr_log2FC_cutoff_info_file_broad=${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_broad.txt
fdr_log2FC_cutoff_info_file_specific=${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_specific.txt

out_num_genes_info_file=${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_num_genes_info.txt
rm -f ${out_num_genes_info_file}

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_broad_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array2
    do
        tissue_of_interest="${array2[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        #filter out genes based on cutoffs
        num_genes=$(awk -F"\t" -v fdr_cutoff="${fdr_cutoff}" -v log2FC_cutoff="${log2FC_cutoff}" '($7!="NA") && ($9!="NA") && ($9<fdr_cutoff) && ($7>log2FC_cutoff || $7< -log2FC_cutoff) {print}' ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt | wc -l)
        echo -e "gtex_v8_${tissue_of_interest}_fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}\t${tissue_of_interest}\tbroad\t${fdr_cutoff}\t${log2FC_cutoff}\t${num_genes}" >> ${out_num_genes_info_file}
    done<${all_broad_tissue_file}
done<${fdr_log2FC_cutoff_info_file_broad}

out_file_path_with_file_header=${gtex_v8_dir}/output/novaseq_12_15_18_denylist_filtered_specific_tissues
while IFS=$'\t' read -a array
do
    fdr_cutoff="${array[0]}"
    log2FC_cutoff="${array[1]}"
    while IFS=$'\t' read -a array2
    do
        tissue_of_interest="${array2[0]}"
        tissue_of_interest="${tissue_of_interest// /_}"

        #filter out genes based on cutoffs
        num_genes=$(awk -F"\t" -v fdr_cutoff="${fdr_cutoff}" -v log2FC_cutoff="${log2FC_cutoff}" '($7!="NA") && ($9!="NA") && ($9<fdr_cutoff) && ($7>log2FC_cutoff || $7< -log2FC_cutoff) {print}' ${out_file_path_with_file_header}_${tissue_of_interest}.protein_coding.with_gene_meta_info.txt | wc -l)
        echo -e "gtex_v8_${tissue_of_interest}_fdr_cutoff_${fdr_cutoff}_log2FC_cutoff_${log2FC_cutoff}\t${tissue_of_interest}\tspecific\t${fdr_cutoff}\t${log2FC_cutoff}\t${num_genes}" >> ${out_num_genes_info_file}
    done<${all_specific_tissue_file}
done<${fdr_log2FC_cutoff_info_file_specific}


#merge with output summary data from gtex
type_iter_header=gtex_v8

grep ${type_iter_header} ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.txt > ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.txt

echo -e "type\tempirical_score_avg\tperm_score_avg_pval\tperm_score_avg_mean\tperm_score_avg_sd\tperm_score_avg_z_score\ttissue_type\tanalysis_type\tfdr_cutoff\tlog2FC_cutoff\tnum_genes" > ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.with_meta_info.txt
join -t $'\t' <(sort -k1,1 ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.txt) -21 <(sort -k1,1 ${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_num_genes_info.txt) >> ${perm_type_results_dir}/all_score_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.with_meta_info.txt

grep ${type_iter_header} ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.txt
grep ${type_iter_header} ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt | awk '{print $NF}' | awk -F"_" '{for(i=1;i<=NF-3;i++){printf $i"_"};printf $(NF-2)"\n"}' > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.ids
grep ${type_iter_header} ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.txt | awk '{print $NF}' | awk -F"_" '{print $NF}' > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.distance

paste <(cat ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.ids) <(cat ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.txt) <(cat ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.distance) > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.revised.txt

echo -e "label\tempirical_count\tcount_pval\tperm_count_mean\tperm_count_sd\tcount_z_score\tclass\ttissue_type\tanalysis_type\tfdr_cutoff\tlog2FC_cutoff\tdistance\tnum_genes" > ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.with_meta_info.txt
join -t $'\t' <(sort -k1,1 ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.revised.txt) -21 <(sort -k1,1 ${gtex_v8_dir}/gtex_v8_fdr_log2FC_cutoff_info_file_num_genes_info.txt) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$9"\t"$14}' >> ${perm_type_results_dir}/all_count_enrichment_hcondel_denylist_filtered_final_stats.${type_iter_header}.with_meta_info.txt

