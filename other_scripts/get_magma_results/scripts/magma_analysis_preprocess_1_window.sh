#preprocess magma files - note that this was part of magma_analysis_v2.sh - got reorganized
#I updated to clean the file - I removed the section "data frames for integrating with original dataset" as that is represented in mpra_set_1_extract_meta_info.sh
#I changed "gene_sig" to "gene_overlap" and "Sig" to "Overlap"
#I removed the LOXL2 analyses - for that go back to magma_analysis_v2.sh 

#*********************** initial parameters ***********************#
filterdir=/cluster_path/ape_project/deletions_project/filtering_new

input_header=magma_analysis

magma_dir=/cluster_path/ape_project/deletions_project/magma
script_dir=${magma_dir}/scripts

#chain file for lifting from hg38 to hg19
liftChainFile=/cluster_path/ape_project/deletions_project/chain_files/hg38ToHg19.over.chain.gz

#store magma program
bin_dir=/cluster_path/bin

#*********************** downloads ***********************#

#download magma program
wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.09a.zip -P ${bin_dir}
unzip ${bin_dir}/magma_v1.09a.zip

mkdir -p ${magma_files}/ref
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip -P ${magma_files}/ref

#European LD
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P ${magma_files}/ref

#*********************** preprocess first ***********************#

#filter out allele frequencies less than 0.01
plink --bfile ${magma_dir}/ref/g1000_eur --maf 0.01 --make-bed --out ${magma_dir}/ref/g1000_eur_maf_filtered
awk '{print $1"\t"$1"_"$4"_"$5"_"$6"\t"$3"\t"$4"\t"$5"\t"$6}' ${magma_dir}/ref/g1000_eur_maf_filtered.bim > ${magma_dir}/ref/g1000_eur_maf_filtered.new.bim
mv ${magma_dir}/ref/g1000_eur_maf_filtered.new.bim ${magma_dir}/ref/g1000_eur_maf_filtered.bim

####### from the ensembl gene info, change chromsome names to ucsc chromsome names ######
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt -P ${magma_dir}/ref

awk -F"\t" '{print $1}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.txt | sort | uniq >  ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.all_chr_names
awk 'NR>1{print}' ${magma_dir}/ref/hg38.chromAlias.txt | awk '{if(($2~/HC/) || ($2~/HG/) || ($2~/HS/)){print $1"\tCHR_"$2}else{print $1"\t"$2}}' | sort -k2,2 > ${magma_dir}/ref/hg38.chromAlias.parsed.1.txt
awk 'NR>1{print}' ${magma_dir}/ref/hg38.chromAlias.txt | awk '{if(($3~/HC/) || ($3~/HG/) || ($3~/HS/)){print $1"\tCHR_"$3}else{print $1"\t"$3}}' | sort -k2,2 > ${magma_dir}/ref/hg38.chromAlias.parsed.2.txt
awk 'NR>1{print}' ${magma_dir}/ref/hg38.chromAlias.txt | awk '{if(($4~/HC/) || ($4~/HG/) || ($4~/HS/)){print $1"\tCHR_"$4}else{print $1"\t"$4}}' | sort -k2,2 > ${magma_dir}/ref/hg38.chromAlias.parsed.3.txt

join -t $'\t' -12 <(cat ${magma_dir}/ref/hg38.chromAlias.parsed.1.txt) -21 <(cat ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.all_chr_names) > ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.1.txt
join -t $'\t' -12 <(cat ${magma_dir}/ref/hg38.chromAlias.parsed.2.txt) -21 <(cat ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.all_chr_names) > ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.2.txt
join -t $'\t' -12 <(cat ${magma_dir}/ref/hg38.chromAlias.parsed.3.txt) -21 <(cat ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.all_chr_names) > ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.3.txt

cat ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.1.txt > ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt
cat ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.2.txt >> ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt
cat ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.3.txt >> ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt

sort -k1,1 ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt > ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.sorted.txt
mv ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.sorted.txt ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt

#sanity check that all chr ids match
#some ids just are not in ${magma_dir}/ref/hg38.chromAlias.txt, see output, thats fine as deletions aren't on those chromosomes anyway
#join -v2 -t $'\t' -11 <(sort -k1,1 ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt) -21 <(cat ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.all_chr_names) 

awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info.txt | awk 'NR>1{print}' | sort -k1,1 -k2,2n > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt

#convert ensembl chr name to match ucsc
 
join -t $'\t' -11 -21 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt >  ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA_temp.txt 

#put ucsc chr name first
paste -d "\t" <(cut -f17  ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA_temp.txt ) <(cut -f2-16 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA_temp.txt ) > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA_revised.txt

#sanity check that all chromosomes are matched with ucsc/ensembl
#join -v1 -t $'\t' -11 -21 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt ${magma_dir}/ref/human_GRCh38_p13_ensembl_hg38_ucsc_id_map.txt

mv ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA_revised.txt ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt 

###### get the TSS positions #########
#chr, TSS, gene name, ensembl id
awk -F"\t" '($4!="NA") {print $1"\t"$4-1"\t"$4"\t"$15"\t"$13}'  ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt | sort -k1,1 -k2,2n | uniq > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_TSS.txt

#sanity check that genes have only one ensembl id
#awk -F"\t" '($1~"[0-9][0-9]*$"){print $13"\t"$15}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt | sort | uniq | awk '{print $1}' | wc -l
#awk -F"\t" '($1~"[0-9][0-9]*$"){print $13"\t"$15}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt | sort | uniq | awk '{print $2}' | wc -l

############### get conservation info as a conditional variable for magma #############
#get gene start and end for all genes
#remove _alt and _fix also as those chr that has funky coordinates, coordinates exceed genome sizes
awk '{if($NF=="1"){print $1"\t"$2-1"\t"$2"\t+\t"$15"\t"$13} else{print $1"\t"$3-1"\t"$3"\t-\t"$15"\t"$13}}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt | sort -k1,1 -k2,2n | uniq | grep -E -v "_alt|_fix" > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord.bed

#calculate prop of regulatory domain having conserved sites
#extend 50 kb upstream to capture conserved regions, use this for conditional analysis later
extendBP_up_regulatory_domain=50000
extendBP_down_regulatory_domain=500
python ${script_dir}/extendCoordinatesV2_upstream_downstream.py ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord.bed ${magma_dir}/ref/hg38.chrom.sizes ${extendBP_up_regulatory_domain} ${extendBP_down_regulatory_domain} ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.bed

awk -F'[\t|]' '{print $4"\t"$5"\t"$6"\t"$4"|"$5"|"$6}'  ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt > ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_human_cons_ids.txt
sort -o -k1,1 -k2,2n ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_human_cons_ids.txt 

#run bedtools coverage to get proportion 
bedtools coverage -a ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.bed -b ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered_with_human_cons_ids.txt > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.cons_coverage.txt

#get gene id and proportion conserved - save for conditional analysis later
awk '{print $9"\t"$10"\t"$11"\t"$13}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.cons_coverage.txt > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.cons_coverage.gene_info.txt

############# get hg19 coordinates because GWAS tend to use hg19 ################

#from list of genes, get genome coordinates for all of them, added 11/15/22 grep -E -v "_alt|_fix"
awk '{print $1"\t"$2-1"\t"$3"\t"$1"|"$2-1"|"$3"|"$15"|"$13"\t0\t"$NF}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_master_info_empty_filled_NA.txt | awk 'BEGIN{FS=OFS="\t"} {gsub("-1","-",$6)}1' | awk 'BEGIN{FS=OFS="\t"} {gsub("1","+",$6)}1' | sort -k1,1 -k2,2n | uniq | grep -E -v "_alt|_fix" > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_coord_map.initial.bed

#now, lift to hg19
#lift human TSS coord to chimp panTro4

liftOver -minMatch=0.001 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_coord_map.initial.bed ${liftChainFile} ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_coord_map.initial.hg19_mapped.txt ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_coord_map.initial.hg19_unmapped.txt

awk -F'[\t|]' '{print $1"\t"$2"\t"$3"\t"$10"\t"$7"\t"$8}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_coord_map.initial.hg19_mapped.txt | sort -k6,6 | uniq > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt 

#added 11/16/22 - did not have extensions of the genes before - had it in the original version
extendBP_up_gene_boundaries=35000
extendBP_down_gene_boundaries=10000
python ${script_dir}/extendCoordinatesV2_upstream_downstream.py ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt ${magma_dir}/ref/hg19.chrom.sizes ${extendBP_up_gene_boundaries} ${extendBP_down_gene_boundaries} ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map_extended_upstream_${extendBP_up_gene_boundaries}_downstream_${extendBP_down_gene_boundaries}.txt

#rearrange coordinates, remove stand and just keep gene name/id
#${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt is the gene location file for magma
awk '{print $9"\t"$1"\t"$2"\t"$3}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map_extended_upstream_${extendBP_up_gene_boundaries}_downstream_${extendBP_down_gene_boundaries}.txt | sort -k1,1 | sed 's#chr##g' > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt

#added 11/14/22 - get the conservation statistics just for the hg19 mapped coordinates, make the GWAS meta file
awk '{print $1}' ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.txt | sort | uniq > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.gene_ids.txt
join -t $'\t' -11 <(sort -k1,1 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.gene_ids.txt) -21 <(sort -k1,1 ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_regulatory_domain}_downstream_${extendBP_down_regulatory_domain}.cons_coverage.gene_info.txt) > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_gene_boundaries}_downstream_${extendBP_down_gene_boundaries}.cons_coverage.gene_info.hg19_mapped.txt

echo -e "gene_id\tcons_num\tcons_bp\tcons_prop" > ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_magma_gene_additional_meta.txt
cat ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding.gene_tss_coord_extended_upstream_${extendBP_up_gene_boundaries}_downstream_${extendBP_down_gene_boundaries}.cons_coverage.gene_info.hg19_mapped.txt >> ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_magma_gene_additional_meta.txt

#******** get the windows meta files ************#
#get the del coordinates
awk 'NR>1{print $11"\t"$12"\t"$13"\t"$1}' ${magma_dir}/ref/novaseq_12_15_18_hCONDEL_cleaned_metatable_submission.txt | sort -k1,1 -k2,2n  > ${magma_dir}/ref/hCONDELs_human_del.txt

######### get genes within window_size base window around deletion #########
#copied and modified from mpra_set_1_extract_meta_info.sh

window_size=50000

bedtools window -w ${window_size} -a ${magma_dir}/ref/hCONDELs_human_del.txt -b ${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_TSS.txt | awk '{print $4"\t"$8"\t"$9}' | sort -k1,1 | uniq  > ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_temp.txt
#put NA in fields that aren't joined 
join -t $'\t' -v 1 -11 <(awk '{print $4"\tNA\tNA"}' ${magma_dir}/ref/hCONDELs_human_del.txt | sort -k1,1) -21 <(sort -k1,1 ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_temp.txt) > ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt

cat ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_temp.txt >> ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt
sort -o -k1,1 ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt  
rm -f ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_temp.txt

#combine gene names that are equidistant from deletion, its possible for a deletion to have multiple window_${window_size} genes
python ${script_dir}/combine_window_bed_entries.py ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_revised.txt
sort -k1,1 ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_revised.txt > ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt
rm -f ${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info_revised.txt

### run magma_analysis_preprocess_v2_window.r to get overlap info for window ###

all_genes_vec_file_path=${magma_dir}/ref/biomart_human_GRCh38_p13_protein_coding_gene_coord_map.gene_ids.txt
hcondel_window_metadata_file_path=${magma_dir}/ref/hCONDELs_human_del_window_${window_size}_gene_info.txt
#size of window used to overlap genes
window_size=50000
additional_info="background_all_genes"
output_file_path_with_header="${magma_dir}/ref/hCONDEL"
write_observed_overlap_genes="Y"
#randomly scramble genes based off of number of genes overlapped, set to "Y" to scramble
write_random_overlap_genes="Y"

#random seed to use
seed_num=5

Rscript ${script_dir}/magma_analysis_preprocess_2_v2_window.r ${all_genes_vec_file_path} ${hcondel_window_metadata_file_path} ${window_size} ${additional_info} ${output_file_path_with_header} ${write_observed_overlap_genes} ${seed_num} ${write_random_overlap_genes}

