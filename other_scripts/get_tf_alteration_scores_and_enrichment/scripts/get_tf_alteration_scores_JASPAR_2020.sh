###### preprocess files ##############

#download rna cell line data from https://www.proteinatlas.org/about/download

unzip -c /cluster_path/ape_project/deletions_project/cell_type_rna/rna_celline.tsv.zip | grep HEK > /cluster_path/ape_project/deletions_project/cell_type_rna/HEK293_rna_celline.tsv

rna_meta_file=/cluster_path/ape_project/deletions_project/cell_type_rna/Encode_RNA_Cell_Type_Table.txt

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
del_dir=/cluster_path/ape_project/deletions_project
out_path=/cluster_path/ape_project/deletions_project/cell_type_rna
af_filter_2=0.00001

scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

#revised 6/25/22, human TSS from ensembl 104, 
awk -F"\t" 'NR>1 {print "chr"$1"\t"$4-1"\t"$4"\t"$13}' ${del_dir}/biomart_human_GRCh38_p13_protein_coding_gene_master_info.txt | sort -k4,4 | uniq > ${out_path}/hg38_all_genes_TSS.txt

#first process encode data only , get in format of gene name, tpm
grep EN ${rna_meta_file} > ${rna_meta_file}.encode

rm -f ${rna_meta_file}.revised

#get expressed genes only
while read -a array
do
	cell_type="${array[0]}"
	file_name="${array[1]}"
	meta_info="${array[2]}"
	meta_info_link="${array[3]}"
	
	#remove version info from gene
	awk 'NR>1{print $1}' ${file_name} | awk -F"." '{print $1}' > ${file_name}.gene_names_revised
	awk 'NR>1{print $6}' ${file_name}  > ${file_name}.tpm
	
	paste -d"\t" ${file_name}.gene_names_revised ${file_name}.tpm > ${file_name}.revised

	echo -e "${cell_type}\t${file_name}.revised\t${meta_info}\t${meta_info_link}" >> ${rna_meta_file}.revised
	
done<${rna_meta_file}.encode

#revise HEK file, put HEK sample in new meta file as well 
awk '{print $1"\t"$5}' /cluster_path/ape_project/deletions_project/cell_type_rna/HEK293_rna_celline.tsv > /cluster_path/ape_project/deletions_project/cell_type_rna/HEK293_rna_celline.tsv.revised

grep HEK ${rna_meta_file} | awk '{print $1"\t"$2".revised\t"$3"\t"$4}' >> ${rna_meta_file}.revised

#revise NPC file, put NPC sample in new meta file as well, added 6/25/22
awk 'NR>1{print $2"\t"$4}' /cluster_path/ape_project/deletions_project/cell_type_rna/hoffman_2017_npc_max_tpm.txt | sort -k1,1 > /cluster_path/ape_project/deletions_project/cell_type_rna/hoffman_2017_npc_max_tpm.txt.revised

grep NPC ${rna_meta_file} | awk '{print $1"\t"$2".revised\t"$3"\t"$4}' >> ${rna_meta_file}.revised

###### filter out files - tpm>1 for expressed ######################

while read -a array
do
	cell_type="${array[0]}"
	file_name="${array[1]}"
	meta_info="${array[2]}"
	meta_info_link="${array[3]}"
	
	#remove version info from gene
	
	awk '$2>=1{print}' ${file_name} | sort -k1,1 > ${file_name}.expressed

	echo ${cell_type}
	wc -l ${file_name}.expressed
	
	if [ "${cell_type}" == "NPC" ];
	then
	    #join by gene name
	    #see notes on why those genes were chosen
	    join -v1 -t $'\t' -11 ${file_name}.expressed -21 <(sort -k1,1 ${out_path}/hoffman_2017_replace_gene_names.txt) > ${file_name}.expressed.ensembl_id_same
	    join -t $'\t' -11 ${file_name}.expressed -21 <(sort -k1,1 ${out_path}/hoffman_2017_replace_gene_names.txt)  | awk '{print $3"\t"$2}' > ${file_name}.expressed.ensembl_id_changed
	    
	    cat ${file_name}.expressed.ensembl_id_same > ${file_name}.expressed
	    cat ${file_name}.expressed.ensembl_id_changed >> ${file_name}.expressed
	    
	    sort -k1,1 ${file_name}.expressed > ${file_name}.expressed.temp
	    mv ${file_name}.expressed.temp ${file_name}.expressed
	    
	    rm -f ${file_name}.expressed.ensembl_id_same 
	    rm -f ${file_name}.expressed.ensembl_id_changed 
	    #put the joined and not joined TFs together and output as final file
    fi
    
	#rearrange file to get chr start pos, end pos, gene name, tpm
	join -t $'\t' -11 ${file_name}.expressed -24 ${out_path}/hg38_all_genes_TSS.txt | awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' | sort -k1,1 -k2,2n | uniq > ${file_name}.expressed.tss_info
	
	wc -l ${file_name}.expressed.tss_info
	
done<${rna_meta_file}.revised


##### merge with human deletion coordinates to get closest expressed gene#############


awk '{print $1"|"$3}' ${post_vcf_filter_dir}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2.txt | awk  -F"|" 'BEGIN{OFS="\t";} {$1=$1}1' > ${post_vcf_filter_dir}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_tab_sep.txt 
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_tab_sep.txt

#get the human only coordinates, the inversion info is for extending the bp
paste <( awk '{print $10"\t"$11"\t"$12"\t"$13}' ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} )   | sort -k1,1 -k2,2n > ${post_vcf_filter_dir}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_human_coord.txt


while read -a array
do
	
	cell_type="${array[0]}"
	file_name="${array[1]}"
	meta_info="${array[2]}"
	meta_info_link="${array[3]}"
	
	#bedtools closest -d -a ${post_vcf_filter_dir}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_human_coord.txt -b ${file_name}.expressed.tss_info | awk '{print }' >  ${out_path}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_human_coord_${cell_type}_closest_expressed_gene.txt
    #added 6/25/22
    #take the max TPM for closest gene if there are ties
    python ${scriptdir}/get_max_tpm.py ${out_path}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_human_coord_${cell_type}_closest_expressed_gene.txt ${out_path}/all_files_af_${af_filter_2}_vcf_checked_all_coord_v2_human_coord_${cell_type}_closest_expressed_gene.max.txt
done<${rna_meta_file}.revised

####### filter out TFs expressed for FIMO analysis later #########################
#for genecard_ID_map_full.txt file, first column is FIMO ID, second column is FIMO TF Name, third column is ENSEMBL TF Name, fourth column is TF Name associated with ENCODE called peaks
meme_suite_files_path=/cluster_path/ape_project/deletions_project/meme_suite_analysis

#how did I get TF map?
#downloaded TF info from biomart
#use GO Term 0003700, https://www.ebi.ac.uk/QuickGO/term/GO:0003700, check uniq to get unique results
#name the file GRCh38.p12_ensembl_TF_map.txt 
#also download all genes without TF filter and name it  ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map.txt  for comparison

awk '{print $1"\t"$NF }' ${meme_suite_files_path}/GRCh38.p12_ensembl_TF_map.txt | sort -k1,1 | uniq >   ${meme_suite_files_path}/GRCh38.p12_ensembl_TF_map_unique.txt
awk '{print $1"\t"$NF }' ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map.txt | sort -k1,1 | uniq >   ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map_unique.txt

##### JASPAR #######
#wget http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt -P ${meme_suite_files_path}
grep MOTIF ${meme_suite_files_path}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt | awk '{print $2"\t"$3}' > ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_1.txt
#replace the "("", and ")"
#http://jaspar.genereg.net/docs/
#The name of the transcription factor. As far as possible, the name is based on the standardized Entrez gene symbols. In the case the model describes a transcription factor hetero-dimer, two names are concatenated, such as RXR-VDR.
awk '{print $2}' ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_1.txt | sed 's#(.*##g' | awk '{print toupper($1)}'> ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_2.txt
paste <(cat ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_1.txt) <(cat ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_2.txt) | awk '{print $1"\t"$2"\t"$3"\t"$3}' > ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt

while read -a array
do

	cell_type="${array[0]}"
	file_name="${array[1]}"
	meta_info="${array[2]}"
	meta_info_link="${array[3]}"
		
	#add gene name to expressed file
	join -t $'\t' -11 ${file_name}.expressed -21 ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map_unique.txt | awk '{print $1"\t"$3}' > ${file_name}.expressed_with_gene_name
	
	#JASPAR
	python ${scriptdir}/filter_expressed_tfs.py ${file_name}.expressed_with_gene_name ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt ${out_path}/${cell_type}_JASPARv2020_expressed_TF.txt
	sort -k1,1 ${out_path}/${cell_type}_JASPARv2020_expressed_TF.txt | uniq > ${out_path}/${cell_type}_JASPARv2020_expressed_TF.temp.txt
	mv ${out_path}/${cell_type}_JASPARv2020_expressed_TF.temp.txt ${out_path}/${cell_type}_JASPARv2020_expressed_TF.txt
	
	#count total number of TFs expressed for info
	join -11 <( cat ${file_name}.expressed | sort -k1,1) -21 ${meme_suite_files_path}/GRCh38.p12_ensembl_TF_map_unique.txt  | awk '{print $NF}' | uniq > ${out_path}/${cell_type}_total_expressed_TF.txt
	
	echo ${cell_type}
	wc -l ${out_path}/${cell_type}_JASPARv2020_expressed_TF.txt
	wc -l ${out_path}/${cell_type}_total_expressed_TF.txt
	
done<${rna_meta_file}.revised


############# run FIMO ##########################################################

#FIMO analysis
#here, I get the macaque coordinate as well, if a macaque coordinate does not exist, then the field will simply be "NA"
#every cell type will have its own file

input_header=cell_type_specific_FIMO_analysis
real_data_out_header=cell_type_specific_FIMO_analysis

del_dir=/cluster_path/ape_project/deletions_project

out_path=/cluster_path_temp/cell_type_specific_rna
scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

species_ref_meta_file=${del_dir}/chimp_ref_human_macaque_input_file.txt

#extend maxBPWindow to capture all surrounding motifs around del site
maxBPWindow=60

#cutOffVal is the pvalue cutoff for fimo
cutOffVal=0.0001


combined_fasta_file=${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id.fa

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering

mkdir -p ${out_path}

#pairwise_comp_file are the pairwise comparisons to do for fimo
pairwise_comp_file=/cluster_path/ape_project/deletions_project/meme_suite_analysis/pairwise_comp_info.all_comp.txt

species_ref_meta_file=${del_dir}/chimp_ref_human_macaque_input_file.txt

#run all cell types
cell_type_file=/cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_cell_list.txt

meme_suite_files_path=/cluster_path/ape_project/deletions_project/meme_suite_analysis


###### TF specific parameters #####

###### JASPARv2020 ##########
memeFileName=${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt
motif_file=/cluster_path/ape_project/deletions_project/meme_suite_analysis/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
tf_db_name=JASPARv2020

################ IMPORTANT CHECK ##################################

#look at encode IDs and see which gene names don't have an ensembl match 

awk '{print $2}' /cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_all_hg19_hg38_fileID_TFname_map.txt | sort | uniq > /cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_all_hg19_hg38_fileID_TFname_map.TFNameOnly.txt

#which encode TFs don't have an ensembl id?
join -v2 -12 -21 <(sort -k2,2 ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map_unique.txt) <(sort -k1,1 /cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_all_hg19_hg38_fileID_TFname_map.TFNameOnly.txt  ) | wc -l
#aside from histone genes, there is only one name that didn't match ZZZ3, which has alias AC118549.1 in ensembl
#AC118549.1 is not in hocomoco though so its ok to be ignored

#which ones in JASPAR don't have a match?
python ${scriptdir}/filter_expressed_tfs.py ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map_unique.txt ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt ${out_path}/JASPARv2020_ensembl_gene_name_check.txt
join -v1 -12 <(sort -k2,2 ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt) -21 <(sort -k1,1 ${out_path}/JASPARv2020_ensembl_gene_name_check.txt) 
#all look like mice genes, EWSR1-FLI1 is a mouse fusion gene it seems https://en.wikipedia.org/wiki/EWS/FLI

#check that all TFs in each database have an ensembl ID
while read -a array
do
	cell_type="${array[0]}"
	file_name="${array[1]}"
	meta_info="${array[2]}"
	meta_info_link="${array[3]}"
	
	rm -f ${file_name}.checked
	
	join -t $'\t' -11 <(sort -k1,1 ${file_name}) -21 <(sort -k1,1 ${meme_suite_files_path}/GRCh38.p12_ensembl_all_genes_map_unique.txt) | sort -k1,1 | awk '{print $1"\t"$3}' | uniq > ${file_name}.ensembl_id.checked
	
	echo -e "${cell_type}"
	echo -e "Number of genes in original gene exp table"
	wc -l ${file_name}
	
	echo -e "Number of genes mapped to my version of ensembl GRCh38.p12"
	wc -l ${file_name}.ensembl_id.checked
	
	#check with each of the TF databases?
    #JASPAR
    echo -e "JASPAR number of TFs not mapped"
    python ${scriptdir}/filter_expressed_tfs.py ${file_name}.ensembl_id.checked ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt ${out_path}/JASPARv2020_ensembl_gene_name_check.exp_table.txt
    join -v1 -12 <(sort -k2,2 ${meme_suite_files_path}/JASPARv2020_genecard_ID_map_full.txt) -21 <(sort -k1,1 ${out_path}/JASPARv2020_ensembl_gene_name_check.exp_table.txt) | wc -l

	#rearrange file to get chr start pos, end pos, gene name, tpm
	#join -t $'\t' -11 ${file_name}.expressed -24 ${out_path}/hg38_all_genes_TSS.txt | awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' | sort -k1,1 -k2,2n | uniq > ${file_name}.expressed.tss_info
done<${rna_meta_file}.revised

###########################################

mpra_set_dir=/cluster_path/ape_project/deletions_project/mpra_1_set

grep -v -E "EMVAR|LCL_HEPG2_POS" ${mpra_set_dir}/MPRA_1_JX_230_BP_36K_oligos.fa | grep "CC" | grep ">" | sed 's/>//g' | awk -F"[|]" '{print $13"\t"$14"\t"$15"\t+\t"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14"|"$15"|"$16"|"$17"|"$18"|"$19"|"$20}' | sed 's/|CC//g'  > ${out_path}/${real_data_out_header}_chimp_coord.txt
grep -v -E "EMVAR|LCL_HEPG2_POS" ${mpra_set_dir}/MPRA_1_JX_230_BP_36K_oligos.fa | grep "HH" | grep ">" | sed 's/>//g' | awk -F"[|]" '{print $16"\t"$17"\t"$18"\t"$19"\t"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14"|"$15"|"$16"|"$17"|"$18"|"$19"|"$20}' | sed 's/|HH//g'  > ${out_path}/${real_data_out_header}_human_coord.txt

#10/23/21 changed to deny list to speed up
#denylist_filtered_ids_file=/cluster_path/ape_project/deletions_project/denylist/novaseq_12_15_18_denylist_filtered_ids.txt

#awk -F"|" '{print $7"\t"$8"\t"$9"\t+\t"$0}' ${denylist_filtered_ids_file} > ${out_path}/${real_data_out_header}_chimp_coord.txt
#awk -F"|" '{print $10"\t"$11"\t"$12"\t"$13"\t"$0}' ${denylist_filtered_ids_file} > ${out_path}/${real_data_out_header}_human_coord.txt

#for additional species other than chimp/human (check species_ref_meta_file), lift to species of interest
#we can also extend bp later

#rm -f ${out_path}/${real_data_out_header}_all_lift_ids.txt

awk '!(($1~ "chimp") || ($1 ~"human")) {print}' ${species_ref_meta_file} > ${out_path}/${real_data_out_header}_additional_species.txt
while read -a array
do
	speciesName=${array[0]}
	
	chrSizeFile=${array[2]}
	liftChainFile=${array[3]}
	liftOver -minMatch=0.001 ${out_path}/${real_data_out_header}_chimp_coord.txt ${liftChainFile} ${out_path}/${real_data_out_header}_${speciesName}_mapped.txt ${out_path}/${real_data_out_header}_${speciesName}_unmapped.txt
	sort -k1,1 -k2,2n ${out_path}/${real_data_out_header}_${speciesName}_mapped.txt > ${out_path}/${real_data_out_header}_${speciesName}_coord.txt 
	rm -f ${out_path}/${real_data_out_header}_${speciesName}_mapped.txt 
	rm -f ${out_path}/${real_data_out_header}_${speciesName}_unmapped.txt
    #awk '{print $5}' ${out_path}/${real_data_out_header}_${speciesName}_coord.txt  >> ${out_path}/${real_data_out_header}_all_lift_ids.txt
    
done<${out_path}/${real_data_out_header}_additional_species.txt



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
	python ${scriptdir}/extendCoordinatesCHIPV3b.py ${out_path}/${real_data_out_header}_${speciesName}_coord.txt ${chrSizeFile} ${maxBPWindow} ${out_path}/${real_data_out_header}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt
        
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
#numCoord=10
#numLines=$(echo ${numCoord} | awk '{print $1*3*2}')
#out_split_path=/cluster_path_temp/roadmap/LCL_data/fimo_coord_split
#mkdir -p ${out_split_path}
#write alternative script to split files
#split -a 3 -l ${numLines} ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa ${out_split_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.

numCoord=10
out_split_path=/cluster_path_temp/cell_type_specific_rna/fimo_coord_split
rm -fr ${out_split_path}
mkdir -p ${out_split_path}
python ${scriptdir}/fimo_split_fasta.py ${out_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa ${numCoord} ${out_split_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.

ls -l ${out_split_path}/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.* | awk '{print $NF}' > ${out_path}/all_fimo_file_names.txt

#run fimo checking for encode annotation

while read -a array
do
	cell_header_id=${array[1]}
	
	#encode file name to overlap
	#0 for overlap
	#1 for no overlap
	#this table is useful eventually for overlapping and looking at actual del sites
	#created from mpradel_encode_mark_build_table_master_script_dist_intensity_combined.sh
	encodeFileName=/cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_${cell_header_id}_closest_dist_count_table_binary.txt
	
	encodeIDMapFileName=/cluster_path/ape_project/deletions_project/tf_misc_files/encode_TF_mpradel_${cell_header_id}_hg19_hg38_fileID_TFname_map.txt
	#TFs expressed in cell type
	#expressedTFsFileName=/cluster_path/ape_project/deletions_project/cell_type_rna/all_species_${cell_header_id}_HOCOMOCOv11_expressed_TF.txt
	expressedTFsFileName=/cluster_path/ape_project/deletions_project/cell_type_rna/all_species_${cell_header_id}_${tf_db_name}_expressed_TF.txt
	
	#humans, chimps, macaques, assume they all express the same TFs from human data, only need to run once per database
	#human
	echo -e "/cluster_path/ape_project/deletions_project/cell_type_rna/${cell_header_id}_${tf_db_name}_expressed_TF.txt" > ${expressedTFsFileName}
	#chimp
	echo -e "/cluster_path/ape_project/deletions_project/cell_type_rna/${cell_header_id}_${tf_db_name}_expressed_TF.txt" >> ${expressedTFsFileName}
	#macaque
	echo -e "/cluster_path/ape_project/deletions_project/cell_type_rna/${cell_header_id}_${tf_db_name}_expressed_TF.txt" >> ${expressedTFsFileName}
	
	#run fimo, centered only, only using expressed genes, not checking for encode annotation
	useEncodeAnnotation=0
	lookAtCenteredOnly=1
	lookAtExpressedTFs=1
	#HOCOMOCO old ID
	#out_header=${cell_header_id}_del_no_encode_annotation_expressed_TF
	out_header=${cell_header_id}_del_no_encode_annotation_expressed_TF_${tf_db_name}
	#for no expressed TF
	#out_header=${cell_header_id}_del_no_encode_annotation_no_expressed_TF_${tf_db_name}
	
	#remove old files, put 7/9/21 4:18pm
	rm -f ${out_path}/${out_header}*
	#v3 takes meme file name as input, same with ${scriptdir}/fimo_run_iter_v4.sh
	#v4 takes an encode ID file to map TFs
	#v5 includes genecard_ID_map.txt to correctly match ids with encode
	bash ${scriptdir}/fimo_run_all_submit_v5.sh ${out_header} ${motif_file} ${pairwise_comp_file} ${lookAtCenteredOnly} ${cutOffVal} ${useEncodeAnnotation} ${encodeFileName} ${encodeIDMapFileName} ${lookAtExpressedTFs} ${expressedTFsFileName} ${memeFileName} ${out_path}

	#head ${out_path}/${out_header}_fimo_all_runs.txt
	#head ${out_path}/${out_header}_fimo_all_runs.log
	#look at all logs
	#ls -l ${out_path}/*_fimo_all_runs.log
	

done<${cell_type_file}


#combine files after finished running

while read -a array
do
	cell_header_id=${array[1]}
	
	#out_header=${cell_header_id}_del_encode_annotation_${tf_db_name}
	#file_out_header=encode_annote_${tf_db_name}
	#bash ${scriptdir}/fimo_combine_all.sh ${out_header} ${file_out_header} ${out_path} ${pairwise_comp_file}
	
	out_header=${cell_header_id}_del_no_encode_annotation_expressed_TF_${tf_db_name}
	file_out_header=no_encode_annote_expressed_TF_${tf_db_name}
	bash ${scriptdir}/fimo_combine_all.sh ${out_header} ${file_out_header} ${out_path} ${pairwise_comp_file}
	
	#sort by date to double check correct files
	echo ${cell_header_id}
	#ls -t --full-time  ${out_path}/${out_header}_* | tail
done<${cell_type_file}

#ls -l ${out_path}/*_fimo_cleaned_output_final.txt
#wc -l ${out_path}/*_fimo_cleaned_output_final.txt
rm -fr ${out_path}/final_fimo_files
mkdir -p ${out_path}/final_fimo_files

cp ${out_path}/*_del_no_encode_annotation_expressed_TF_JASPAR*_fimo_cleaned_output_final.txt ${out_path}/final_fimo_files
cp ${out_path}/*_del_encode_annotation_JASPAR*_fimo_cleaned_output_final.txt ${out_path}/final_fimo_files

cp ${out_path}/NPC_del_no_encode_annotation_expressed_TF_*_fimo_cleaned_output_final.txt ${out_path}/final_fimo_files
########## remove files ######################

bash ${scriptdir}/fimo_remove_all.sh ${out_path}

