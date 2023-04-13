#this script is the followup from analyze_vcf_stats_post_call_new.sh 
#it filters out coordinates by looking at the underlying MAF alignment
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1
mkdir -p ${post_vcf_filter_dir}

chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3.txt

#combine all files into one master file
#from the deletions that pass the initial filter, screen against rest
rm -f ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion.txt
rm -f ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved.txt
rm -f ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord.txt
rm -f ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord.txt

#we reformat all files so that the contents are now: 
#ORIGNAL conserved sequence, with its human region coord,  ORIGINAL deletion coordinate, gap deletion coord, invInfo 
sgdp_output_dir=/cluster_path_temp/sgdp
while read -a array 
do
    genome_file_id=${array[0]}
    vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
    #remove the hybrid coordinates (both the conserved and deleted coordinates)
    cut -f7-19 ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt | sort | uniq >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion.txt 
    #can't have multiple deletions overlapping same conserved site, so no need to uniq below
    #change type to put  conserved sequence first (along with its human coord), then deleted sequence second, end is the inversion info
    awk '{print $13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$19}' ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved.txt  

    #### for partials ########
 
    #get the union coordinates  (remember we unioned partial deletion coordinates and conserved coordinates)
    awk '{print $7"|"$8"|"$9}'  ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt | sort | uniq  >>  ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord.txt

    #get the deletions
    cut -f13-19  ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt | sort | uniq  >>  ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord.txt 
    
done<${chimp_ref_file_names}

#sort files

sort -k1,1 ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord.txt > ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord_sorted.txt

mv ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord_sorted.txt ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord.txt

#sort -k1,1 -k2,2n ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord.txt > ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord_sorted.txt

#mv ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord_sorted.txt ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord.txt

#${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed is from the finalFilteringProtocol_v4.sh , and contains all the partial conserved coordinates without having done any filtering yet
#remember that file format is: ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed
#union chimp coordinates, union human coordinates,  orig del coordinates, human coord, inv info, chimp conserved coord, human conserved coord

awk '{print $1"|"$2"|"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}' ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed | sort | uniq > ${post_vcf_filter_dir}/chimp_partial_deletion_coordinates_cleaned_union_coord_id_orig_cons_coord.txt 

#from the union coordinates, extract the original conserved coordinates
join -11 ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_union_cons_uniq_coord.txt  -21  ${post_vcf_filter_dir}/chimp_partial_deletion_coordinates_cleaned_union_coord_id_orig_cons_coord.txt | awk 'BEGIN{OFS="\t";} {$1=$1}1' | cut -f2- | uniq  >  ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_cons_uniq_coord.txt

#need to redo intersect of partial conserved region with human specific deletions
bedtools intersect -a ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_cons_uniq_coord.txt -b ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_pass_align_del_uniq_coord.txt -loj > ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion.txt 
 

#get only the conserved coordinates from the partial deletions
cut -f1-3 ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion.txt | sort | uniq > ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion_cons_coord_only.txt
cut -f1-3 ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved.txt | sort | uniq > ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved_cons_coord_only.txt
cut -f1-3 ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion.txt | sort | uniq > ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_cons_coord_only.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion_cons_coord_only.txt > ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_coord_only.txt 
awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved_cons_coord_only.txt  >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_coord_only.txt
awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion_cons_coord_only.txt >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_coord_only.txt

#after combining all coordinates together, get the relevant MAFs in each region
#for conserved region overlapping deletion, get the overarching conserved region
#for deletions overlapping conserved regions, get only the conserved region
#for partial deletions overlapping conserved regions, get the conserved region

 
#input file should be list of all conserved sequences 

mafdir=${sgdp_output_dir}/all_maf_files
rm -fr ${mafdir}
mkdir -p ${mafdir}
refdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way
mafsInRegion -outDir ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_coord_only.txt ${mafdir} ${refdir}/target_panTro4_multiz.maf


#check all coordinates at once 
#windowSize is only relevant for conserved sequences overlapping deletions

windowSize=5

#this is compariing the deletion sequences with the MAF sequence
#file input is conserved region, then gap location,

outputFile=${post_vcf_filter_dir}/all_af_${af_filter}_maf_check_coord_output.txt

awk '{print $0"\t""C"}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_conserved_overlapping_deletion.txt > ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt
awk '{print $0"\t""D"}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_deletion_overlapping_conserved.txt  >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt
awk '{print $0"\t""P"}'  ${post_vcf_filter_dir}/all_files_af_${af_filter}_partial_conserved_overlapping_deletion.txt >> ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt

#final file input is the ORIGNAL conserved sequence, with its human region coord,  ORIGINAL gap coordinate, gap human coord, invInfo
#if deletion is greater than cons coord or goes past it, then mark it down
#I should be able to pull directly from this file to design oligos
python checkDeletionOutputV3.py ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt ${mafdir} ${outputFile} ${windowSize}

num_fields=$(awk '{print NF}' ${post_vcf_filter_dir}/all_files_af_${af_filter}_maf_check_coord.txt | head -1)

#master filter

#no filter
awk '{print}' ${outputFile} | cut -f1-${num_fields} > ${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt

#interesting sites filter
#check that conserved site is present in at least 3 primates + that the overall amount of bases that align to the exact deletion position /total bases in region==0.9
awk '($43>=3)&& ($37==0) &&  ($35/($38*$43))>0.8 {print}' ${outputFile} | cut -f1-${num_fields} > ${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_primate_filter.txt

#get the subset which has macaque sequences
awk '($43>=3)&& ($37==0) &&  ($35/($38*$43))>0.8 && ($15/($18*$23)==1) {print}' ${outputFile} | cut -f1-${num_fields}  >  ${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_macaque_filter.txt

######################

#awk '($43>=3)&& ($35/($(NF-2)*$43))>0.8 {print}' ${outputFile} | cut -f1-${num_fields}  >  ${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt
awk '($43>=3)&& ($35/($(NF-2)*$43))>0.8 && ($15/($(NF-2)*$22)==1) {print}' ${outputFile} | cut -f1-${num_fields}  >  ${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered_macaque_filter.txt
<<COMMENT
#some stats
#exact base match in window, changing number of required primates
#39 is window, 35 is exact deletion position
awk '($43>=3)&& ($39/($42*$43))==1 {print }' ${outputFile} | wc -l

awk '($43>=4)&& ($39/($42*$43))==1 {print }' ${outputFile} | wc -l

awk '($43==5)&& ($39/($42*$43))==1 {print }' ${outputFile} | wc -l

awk '($43>=3)&& ($39/($42*$43))>=0.8 {print }' ${outputFile} | wc -l

awk '($43>=4)&& ($39/($42*$43))>=0.8 {print }' ${outputFile} | wc -l

awk '($43==5)&& ($39/($42*$43))>=0.8 {print }' ${outputFile} | wc -l



#
awk '($43>=3)&& ($35/($38*$43))==1 {print }' ${outputFile} | wc -l

awk '($43>=4)&& ($35/($38*$43))==1 {print }' ${outputFile} | wc -l

awk '($43==5)&& ($35/($38*$43))==1 {print }' ${outputFile} | wc -l

awk '($43>=3)&& ($35/($38*$43))>=0.8 {print }' ${outputFile} | wc -l

awk '($43>=4)&& ($35/($38*$43))>=0.8 {print }' ${outputFile} | wc -l

awk '($43==5)&& ($35/($38*$43))>=0.8 {print }' ${outputFile} | wc -l


#distribution of sizes filtering for window
awk '($43>=3)&& ($39/($39*$43))>=0.8 {print $3-$2}' ${outputFile} | sort -k1,1n | uniq -c

awk '($43>=3)&& ($39/($39*$43))>=0.8 {print $9-$8}' ${outputFile} | sort -k1,1n | uniq -c

#distribution of sizes filtering for exact deletion position
awk '($43>=3)&& ($35/($38*$43))>=0.8 {print $3-$2}' ${outputFile} | sort -k1,1n | uniq -c

awk '($43>=3)&& ($35/($38*$43))>=0.8 {print $9-$8}' ${outputFile} | sort -k1,1n | uniq -c

#ensuring that macaque also has the base
awk '($43>=3)&& ($35/($38*$43))>0.8 && ($15/($38*$22)==1) {print}' ${outputFile}  | wc -l
COMMENT
