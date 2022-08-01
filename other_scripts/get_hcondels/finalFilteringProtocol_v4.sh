######################### new filtering protocol ###################################
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
mkdir -p ${filterdir}
chimp_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp
coord_dir=/cluster_path/ape_project/deletions_project/deletion_coordinates
#gap filter is the proportion of the sequence: conserved or deletion that contains gaps, we remove any with significant number of gaps
gapFilter=0.01

#bedtools intersect -wo -a ${coord_dir}/initial_net_hdels.bed -b ${coord_dir}/mostConserved_meta.bed | grep "chr1|88387999|88388194|195|chr1|87604124|0|+"
#grep "chr1|88388087|88388164" ${filterdir}/mostConserved_meta_gap_remove_ids.bed

#### clean for gaps from deletions, remove gaps which comprise more than ${gapFilter} of deletion sequence ####
bedtools intersect -wo -a ${coord_dir}/initial_net_hdels.bed -b ${chimp_ref_dir}/gap_coordinates.txt > ${filterdir}/chimp_deletion_gap_overlap.txt
awk '{sumgaps[$1"\t"$2"\t"$3"\t"$4] += $NF;sumlen[$1"\t"$2"\t"$3"\t"$4]=$3-$2}; END{ for (id in sumgaps) { print id"\t"sumlen[id]"\t"sumgaps[id]"\t"sumgaps[id]/sumlen[id] } }' ${filterdir}/chimp_deletion_gap_overlap.txt | awk -v gapFilter=${gapFilter} '($NF>=gapFilter){print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq > ${filterdir}/initial_net_hdels_gap_remove_ids.bed
join -v1 -14 -24 -o'1.1 1.2 1.3 1.4' <(sort -k4,4 ${coord_dir}/initial_net_hdels.bed ) <(sort -k4,4 ${filterdir}/initial_net_hdels_gap_remove_ids.bed) | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${filterdir}/initial_net_hdels_gap_cleaned.bed

#clean for gaps from conserved sequences, remove gaps which comprise more than ${gapFilter} of conserved sequence
#UPDATE: I didn't actually filter out any gaps from conserved sequences since lifOverConservedTest.sh ignores  ${filterdir}/mostConserved_meta_gap_cleaned.bed LOL
bedtools intersect -wo -a ${coord_dir}/mostConserved_meta.bed  -b ${chimp_ref_dir}/gap_coordinates.txt > ${filterdir}/mostConserved_meta_gap_overlap.txt
awk '{sumgaps[$1"\t"$2"\t"$3"\t"$4] += $NF;sumlen[$1"\t"$2"\t"$3"\t"$4]=$3-$2}; END{ for (id in sumgaps) { print id"\t"sumlen[id]"\t"sumgaps[id]"\t"sumgaps[id]/sumlen[id] } }' ${filterdir}/mostConserved_meta_gap_overlap.txt | awk -v gapFilter=${gapFilter} '($NF>=gapFilter){print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq  > ${filterdir}/mostConserved_meta_gap_remove_ids.bed
join -v1 -14 -24 -o'1.1 1.2 1.3 1.4' <(sort -k4,4 ${coord_dir}/mostConserved_meta.bed ) <(sort -k4,4 ${filterdir}/mostConserved_meta_gap_remove_ids.bed) | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${filterdir}/mostConserved_meta_gap_cleaned.bed

#testing
#bedtools intersect -wo -a ${filterdir}/mostConserved_meta_gap_cleaned.bed -b ${chimp_ref_dir}/gap_coordinates.txt > ${filterdir}/test.txt
#awk '{sumgaps[$1"\t"$2"\t"$3"\t"$4] += $NF;sumlen[$1"\t"$2"\t"$3"\t"$4]=$3-$2}; END{ for (id in sumgaps) { print id"\t"sumlen[id]"\t"sumgaps[id]"\t"sumgaps[id]/sumlen[id] } }' ${filterdir}/test.txt | awk -v gapFilter=${gapFilter} '($NF>=gapFilter){print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq  | head

#liftover conserved sites, from panTro4 to hg38
#the 0.001 setting ensures that at least conserved sequence * 0.001 bp need to be mapped for it to not be considered deleted
#cp ../filtering/liftOverConservedTest.sh .
qsub ${filterdir}/liftOverConservedTest.sh 0.001 ${filterdir}/mostConserved_meta_gap_cleaned.bed

#intersection of conserved mapped sites to human reference to see which ones map to same place - these are repetitive regions, don't keep these, and use for swapping
#this also puts the chimp region first (first 3 coordinates), human region afterwards (next 3 coordinates)
bedtools intersect -loj -a ${filterdir}/mostConservedMapped_new_0.001.txt -b ${filterdir}/mostConservedMapped_new_0.001.txt |awk '{print $4}' | awk  -F '|' '{print $1"\t"$2"\t"$3}' | sort | uniq -c | awk '$1==1{print $2"|"$3"|"$4}' > ${filterdir}/mostConservedMapped_new_unique_ids_0.001.txt 

join -14 -21 <(sort -k4,4 ${filterdir}/mostConservedMapped_new_0.001.txt) <(sort ${filterdir}/mostConservedMapped_new_unique_ids_0.001.txt ) | awk -F'[| ]'  '{print $1"\t"$2"\t"$3"\t"$4"|"$5"|"$6}' > ${filterdir}/mostConservedMapped_new_unique_0.001_regions.txt


###########################################
#for both below, I use a 5% cutoff
#downstream, I only lose a few sites setting this cutoff, I use this cutoff to ensure that the human sequence is not substantially
#larger than the deleted sequence

# for mapped conserved regions remove deletions or conserved sites where the human region is much larger than the deleted region  
#this indicates insertion 
awk -F'[\t|]' '($6-$5)/($3-$2)<1.05 {print}' ${filterdir}/mostConservedMapped_new_unique_0.001_regions.txt >  ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt 
#if we want no filtering, simply move files
#mv ${filterdir}/mostConservedMapped_new_unique_0.001_regions.txt ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt 

#for deleted regions, remove regions where human sequence is much larger than deleted sequence 

awk -F'[\t|]' '$(NF-1)/($3-$2)<1.05 {print}' ${filterdir}/initial_net_hdels_gap_cleaned.bed > ${filterdir}/initial_net_hdels_gap_cleaned_human_size_filtered.bed
#if we want no filtering, simply move files
#mv ${filterdir}/initial_net_hdels_gap_cleaned.bed ${filterdir}/initial_net_hdels_gap_cleaned_human_size_filtered.bed
###########################################################################################################


#get the completely deleted regions#, the ones marked as #Deleted in new, as well as the partial deletions
#wc -l ${filterdir}/mostConservedDeleted_new_0.001.txt 3092
awk 'BEGIN { RS = "#"; FS="\n";  ORS="\n"} ; !($1 ~/Duplicated/) {print $2}' ${filterdir}/mostConservedNotMapped_new_0.001.txt | awk 'NR>1{print}' > ${filterdir}/mostConservedDeleted_new_0.001.txt 

#combine all the conserved sequences (deleted or mapped to hg38) into a single file
awk '{print $0"|""Deleted"}' ${filterdir}/mostConservedDeleted_new_0.001.txt > ${filterdir}/allConserved_filter_0.001.txt
cat ${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt  >> ${filterdir}/allConserved_filter_0.001.txt

#with this final set of conserved sites just created, check to see which ones overlap deletion sites 
bedtools intersect -loj -a ${filterdir}/initial_net_hdels_gap_cleaned_human_size_filtered.bed  -b ${filterdir}/allConserved_filter_0.001.txt | awk '$6!=-1{print}' > ${filterdir}/chimp_deletion_mostconserved_coordinates.bed

#remove deletions which "remove" the same human sequence (these deletions map to multiple human sequences, if we keep these, we don't know which sequence to swap in human swap) - these are repetitive sequences in chimp which removed a corresponding human sequence?
#uniq in first line since a deletion can overlap multiple conserved sites
#I only see two such deletions
awk -F'[\t|]' '{print $8"\t"$9"\t"$9+$10"\t"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11}' ${filterdir}/chimp_deletion_mostconserved_coordinates.bed  | sort | uniq >  ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only_coord.bed 
#intersect the human coordinates that the conserved sequences map to with each other, remove sites which appear multiple times
#this indicates that multiple deletions remove the same human sequence
bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only_coord.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only_coord.bed |  awk '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n | uniq -c  |  awk '$1>1{print $5}' > ${filterdir}/deletions_removing_same_humans_seq_ids.txt
join -v1 -14 -21 -o'1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8' <(sort -k4,4 ${filterdir}/chimp_deletion_mostconserved_coordinates.bed) <(sort ${filterdir}/deletions_removing_same_humans_seq_ids.txt) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned.bed

#next revise coordinates for deletion site partially overlapping conserved${filterdir}/mostConservedMapped_new_0.001_human_size_filtered.txt  site
#awk '!(((($2>$6) && ($3<$7)) ||  (($2==$6) && ($3<$7)) || (($2>$6) && ($3==$7))) || (($2<=$6) && ($3>=$7))) {print}' ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq > ${filterdir}/chimp_partial_deletion_coordinates.bed
awk '!(((($2>$6) && ($3<$7)) ||  (($2==$6) && ($3<$7)) || (($2>$6) && ($3==$7))) || (($2<=$6) && ($3>=$7)))  {print}' ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned.bed | awk '{print}' > ${filterdir}/chimp_partial_deletion_coordinates_all_info.bed
del_dir=/cluster_path/ape_project/deletions_project

######### from here below, this is where finalFilteringProtocol v4 differs from v3 ##########################

awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$9+$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' ${filterdir}/chimp_partial_deletion_coordinates_all_info.bed > ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed
#we union partial deletion sites and conserved sites into one coordinate
#output is union chimp coordinates, union human coordinates,  orig del coordinates, human coord, inv info, chimp conserved coord, human conserved coord
python ${del_dir}/getRevisedConservedCoordFromPartialUnion.py ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed

#get ids of partial deletions and conserved sequences for removal from the other files
awk '{print $1"|"$2"|"$3}' ${filterdir}/chimp_partial_deletion_coordinates_all_info.bed | sort | uniq > ${filterdir}/chimp_partial_deletion_coordinates_all_info_partial_only_deletion_ids.txt
awk '{print $5"|"$6"|"$7}' ${filterdir}/chimp_partial_deletion_coordinates_all_info.bed | sort | uniq > ${filterdir}/chimp_partial_deletion_coordinates_all_info_partial_only_conserved_ids.txt


######## prepare for creating final files ##############
#extract the original conserved and deleted elements

awk '{print $1"\t"$2"\t"$3"\t"$4}'  ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned.bed | sort | uniq | awk -F'[\t|]' '{print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$9+$10"\t"$11"\t"$1"|"$2"|"$3}'  > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed
#file output:  chimp deletion region, corresponding human region, invInfo, id for matching later
awk '{print $5"\t"$6"\t"$7"\t"$8}'  ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned.bed | sort | uniq | awk -F'[\t|]' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"|"$2"|"$3}' > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only.bed
#file output:  chimp conserved region, corresponding human region, id for matching later

#first 6 coordinates contain union coordinates on chimp, along with corresponding human coordinates
cut -f1-6 ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed | sort | uniq > ${filterdir}/chimp_partial_deletion_coordinates_cleaned_coord_only.bed

#remove coordinates that are part of partial deletions, these areas will be treated seperately
join -v1 -18 -21 -o'1.1 1.2 1.3 1.4 1.5 1.6 1.7' <(sort -k8,8 ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed ) <(sort -k1,1 ${filterdir}/chimp_partial_deletion_coordinates_all_info_partial_only_deletion_ids.txt | uniq) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only_no_partial.bed

join -v1 -17 -21 -o'1.1 1.2 1.3 1.4 1.5 1.6' <(sort -k7,7 ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only.bed ) <(sort -k1,1 ${filterdir}/chimp_partial_deletion_coordinates_all_info_partial_only_conserved_ids.txt | uniq) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed

#intersect the conserved coordinates with deletion coordinates
bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only_no_partial.bed > ${filterdir}/conserved_overlapping_deletion.txt
#output conserved coord, human coord, deletion coord, human coord, inv info

#intersect the deletion coordinates with conserved coordinates
bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only_no_partial.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$7}' > ${filterdir}/deletion_overlapping_conserved.txt
# deletion coord, human coord, conserved coord, human coord, inv info

#intersect the union deletion and chimp conserved coordinates with deletion coordinates
bedtools intersect -loj -a ${filterdir}/chimp_partial_deletion_coordinates_cleaned_coord_only.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed > ${filterdir}/partial_conserved_overlapping_deletion.txt 

#extract the sites of where conserved sites overlap deletion sites
#note, I don't include the scenario  $2==$8 && $3==$9, where deletion and conserved sites are equal, this is captured in the deletion overlapping conserved file
awk '((($2<$8) && ($3>$9)) ||  (($2==$8) && ($3>$9)) || (($2<$8) && ($3==$9)))  {print }' ${filterdir}/conserved_overlapping_deletion.txt > ${filterdir}/cleaned_conserved_overlapping_deletion.txt

#extract deletion site completely overlapping conserved site
awk '($8!=-1) && ($2<=$8) && ($3>=$9) {print}' ${filterdir}/deletion_overlapping_conserved.txt > ${filterdir}/cleaned_deletion_overlapping_conserved.txt

#wc -l ${filterdir}/partial_conserved_overlapping_deletion.txt
#wc -l ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt 
cut -f1-13 ${filterdir}/partial_conserved_overlapping_deletion.txt > ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt 

#combine all coordinates into one file and get unique coordinates, then sort based on human coordinates
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord.txt
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort | uniq >> ${filterdir}/final_chimp_human_hybrid_replace_coord.txt
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt | sort | uniq >> ${filterdir}/final_chimp_human_hybrid_replace_coord.txt

#remember that the first six coordinates here, for partials, is the union of partial deletion + conserved coordinate
sort -k 4,4 -k 5,5n -k 6,6n ${filterdir}/final_chimp_human_hybrid_replace_coord.txt > ${filterdir}/final_chimp_human_hybrid_replace_coord_sorted.txt
mv ${filterdir}/final_chimp_human_hybrid_replace_coord_sorted.txt ${filterdir}/final_chimp_human_hybrid_replace_coord.txt

#check that all sites only have one strand
#awk '{print $1"\t"$2"\t"$3}' ${filterdir}/final_chimp_human_hybrid_replace_coord.txt | sort | uniq | wc -l
#for partials, below code shouldn't be necessary? - sanity check
#awk '((($2<$8) && ($3>$9)) ||  (($2==$8) && ($3>$9)) || (($2<$8) && ($3==$9)))  {print }' ${filterdir}/partial_conserved_overlapping_deletion.txt  > ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt 

############## stats ############

#deletion overlapping conserved size range
#awk '{print $3-$2"\t"$9-$8 }' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1n

#number of conserved sites total
awk '{print $1"\t"$2"\t"$3}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1n | uniq | wc -l

#number of deletion sites total
awk '{print $1"\t"$2"\t"$3}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1n | uniq | wc -l
####################
#kind of keep old files for reference 
mv ${filterdir}/final_chimp_human_hybrid_replace_coord.txt ${filterdir}/final_chimp_human_hybrid_replace_coord.old.txt
mv ${filterdir}/cleaned_conserved_overlapping_deletion.txt ${filterdir}/cleaned_conserved_overlapping_deletion.old.txt
mv ${filterdir}/cleaned_deletion_overlapping_conserved.txt ${filterdir}/cleaned_deletion_overlapping_conserved.old.txt
##################### liftOver check if human coord matches right coordinates############################

#in order to check that the human coordinates are annotated correctly, we go through each of the checks

########### for conserved sequences ######################

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 

############## for deletion coordinates ####################

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2-1"\t"$3+1"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2+1"\t"$3-1"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 
#################### for partial coordinates ###########################

#for deletion sequence lying right of conserved sequence
python ${del_dir}/getRevisedConservedCoordFromPartialUnionCheck.py ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed ${filterdir}/chimp_partial_deletion_coordinates_test.bed
awk '$7=="R"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14}' ${filterdir}/chimp_partial_deletion_coordinates_test.bed > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2"\t"$3+1"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3-1"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 


awk '$7=="Rinv"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14}' ${filterdir}/chimp_partial_deletion_coordinates_test.bed > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2"\t"$3+1"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2+1"\t"$3"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 

###for deletion sequence lying to left of conserved sequence###

python ${del_dir}/getRevisedConservedCoordFromPartialUnionCheck.py ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed ${filterdir}/chimp_partial_deletion_coordinates_test.bed
awk '$7=="L"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14}' ${filterdir}/chimp_partial_deletion_coordinates_test.bed > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2-1"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2+1"\t"$3"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 



awk '$7=="Linv"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14}' ${filterdir}/chimp_partial_deletion_coordinates_test.bed > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2-1"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3-1"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 

#### for deletion with conserved sequence on both sides #####


python ${del_dir}/getRevisedConservedCoordFromPartialUnionCheck.py ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed ${filterdir}/chimp_partial_deletion_coordinates_test.bed
awk '$7=="0"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14}' ${filterdir}/chimp_partial_deletion_coordinates_test.bed > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 


######## all partials ##########

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}'   ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt > ${filterdir}/check_file.txt
testfile=${filterdir}/check_file.txt
${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt 

chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${testfile} | sort -k1,1 -k2,2n | uniq >  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt
liftOver -minMatch=0 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_unmapped.txt

awk -F"[\t|]" '{print}' ${testfile} | sort -k1,1 -k2,2n | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$8}' ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped.txt  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 

diff ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed  ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed  
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 


###########
###### followup liftover tests from above, below is an example where extension of even a single bp of conserved region could find some deletions
grep 68632694 ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed
grep 68632694 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed 
grep 68632694 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 
liftOver -minMatch=0.001 ${del_dir}/lifttest.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${del_dir}/lifttest_mapped.txt ${del_dir}/lifttest_unmapped.txt
head ${del_dir}/lifttest_mapped.txt 

grep 7667246 ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed
grep 7667246 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed 
grep 7667246 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 
liftOver -minMatch=0.001 ${del_dir}/lifttest2.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${del_dir}/lifttest2_mapped.txt ${del_dir}/lifttest2_unmapped.txt
head ${del_dir}/lifttest2_mapped.txt


grep 92726196 ${filterdir}/chimp_partial_deletion_coordinates_all_info_parsed.bed
grep 92726196 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_orig.bed 
grep 92726196 ${filterdir}/final_chimp_human_hybrid_replace_coord_lift_test_mapped_check_lift.bed 
liftOver -minMatch=0 ${del_dir}/lifttest3.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${del_dir}/lifttest3_mapped.txt ${del_dir}/lifttest3_unmapped.txt
head ${del_dir}/lifttest3_mapped.txt


######################################################################
############# liftOverCheck #######################
chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38

#why is it partially deleted from new and not completely deleted? is chain file different from net?
#use ${filterdir}/chimp_deletion_mostconserved_coordinates.bed to check if the strand is same as liftover

awk -F"[\t|]" '{print $1"\t"$2-1"\t"$3+1"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${filterdir}/chimp_deletion_mostconserved_coordinates.bed | sort -k1,1 -k2,2n | uniq >  ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test.bed
liftOver -minMatch=0 ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test.bed  ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped.bed ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_unmapped.bed

awk -F"[\t|]" '{print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$9+$10"\t"$11}' ${filterdir}/chimp_deletion_mostconserved_coordinates.bed  | sort -k1,1 -k2,2n | uniq  > ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped_diff_check_orig.bed 
awk -F"[\t|]" '{print $4"\t"$5"\t"$6"\t"$1"\t"$2+1"\t"$3-1"\t"$8}' ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped.bed | sort -k1,1 -k2,2n | uniq > ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped_diff_check_lift.bed

diff ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped_diff_check_orig.bed  ${filterdir}/chimp_deletion_mostconserved_coordinates_lift_test_mapped_diff_check_lift.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_orig.bed 
rm -f ${filterdir}/chimp_partial_deletion_coordinates_cleaned_lift_test_mapped_diff_check_lift.bed 

#take union of deletion and conserved sequence

#get ids of partial coordinates, then use that to completely remove conserved coordinate with partial deletions
#########################################################################################################

#put in the new partial deleted revised conserved coordinates with the conserved sequences that did not partially overlap deletion sites
#awk -F'[\t|]' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised.bed > ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised_partial_only.bed 
#partial cons coord, original conserved coord, human mapping coord


#REDO the intersect, intersect deleted sites with conserved regions, conserved regions with deleted regions, 
#this time the conserved sequences contain only the parts that overlap completely with the deleted sequences
#also rearrange columns

#intersect conserved sequences with deleted sequences

#bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed  > ${filterdir}/conserved_overlapping_deletion.txt
#bedtools intersect -loj -a ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised_partial_only.bed  -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6}' >>  ${filterdir}/conserved_overlapping_deletion.txt
#conserved coord, human coordinate, deleted region - primary columns
#human coord corresponding to deleted region, original conserved coord (if partial) - secondary columns

#and also deletion regions with conserved sites
#bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5+$6"\t"$7"\t"$8"\t"$9"\t.\t.\t.\t"$10"\t"$11"\t"$12}' > ${filterdir}/deletion_overlapping_conserved.txt
#bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed -b ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised_partial_only.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5+$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' >> ${filterdir}/deletion_overlapping_conserved.txt
#deletion coord, human coord, conserved region - primary columns
#original conserved coord, human coord - secondary columns

#secondary columns differ between the two files

#extract the sites of where conserved sites overlap deletion sites
#note, I don't include the scenario  $2==$8 && $3==$9, where deletion and conserved sites are equal, this is captured in the deletion overlapping conserved file
awk '((($2<$8) && ($3>$9)) ||  (($2==$8) && ($3>$9)) || (($2<$8) && ($3==$9)))  {print }' ${filterdir}/conserved_overlapping_deletion.txt > ${filterdir}/cleaned_conserved_overlapping_deletion.txt

#extract deletion site completely overlapping conserved site
#get the deletions partially overlapping conserved regions
awk '($8!=-1) && ($2<=$8) && ($3>=$9) {print}' ${filterdir}/deletion_overlapping_conserved.txt > ${filterdir}/cleaned_deletion_overlapping_conserved.txt


#combine all coordinates into one file and get unique coordinates, then sort based on human coordinates
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord.txt
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort | uniq >> ${filterdir}/final_chimp_human_hybrid_replace_coord.txt

sort -k 4,4 -k 5,5n -k 6,6n ${filterdir}/final_chimp_human_hybrid_replace_coord.txt > ${filterdir}/final_chimp_human_hybrid_replace_coord_sorted.txt
mv ${filterdir}/final_chimp_human_hybrid_replace_coord_sorted.txt ${filterdir}/final_chimp_human_hybrid_replace_coord.txt

################################################################

#get information on whether or not sequence is inverted, this will be a meta file used for other programs where inversion info is important
#I cheat by using liftover again to get this position.. too much of a hassle to go back to beginning of code
#remember that for conserved sequence being completely deleted, we have no idea what the orientation is
#so its encoded as 0 in the meta file outputted below

filterdir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38
awk '{print $1"\t"$2-1"\t"$3+1"\t"$1"|"$2"|"$3"\t""0""\t""+"}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort | uniq > ${filterdir}/cleaned_deletion_overlapping_conserved_temp.txt
liftOver -minMatch=0.001 ${filterdir}/cleaned_deletion_overlapping_conserved_temp.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/cleaned_deletion_overlapping_conserved_temp_mapped.txt ${filterdir}/cleaned_deletion_overlapping_conserved_temp_unmapped.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3"\t""0""\t""+" }' ${filterdir}/final_chimp_human_hybrid_replace_coord.txt | uniq > ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp.txt
liftOver -minMatch=0.001 ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp_mapped.txt ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp_unmapped.txt
#+ will be positive strand, - will be minus strand, 0 will be unknown (since entire sequence is deleted)
awk '{print $4"\t"$NF}' ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp_mapped.txt > ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_inv_info.txt
awk '!($0~"#"){print $4"\t"0}' ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_temp_unmapped.txt >> ${filterdir}/final_chimp_human_hybrid_replace_coord_conserved_inv_info.txt

#final check if any of the human coordinates intersect with each other
#awk '{print $4"\t"$5"\t"$6}' ${filterdir}/final_chimp_human_hybrid_replace_coord.txt > ${filterdir}/final_human_only_coord_check.txt
#bedtools intersect -loj -a ${filterdir}/final_human_only_coord_check.txt -b ${filterdir}/final_human_only_coord_check.txt | awk '(($2!=$5) || ($3!=$6)) {print}' | awk '($2!=$5) {print}' | awk '($3!=$6) {print}' | head
#############################################

#get the fasta file
chimp_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp 

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${filterdir}/final_chimp_human_hybrid_replace_coord.txt > ${filterdir}/final_chimp_human_hybrid_replace_coord_get_seq.txt
bedtools getfasta -name -fi ${chimp_ref_dir}/panTro4.fa -bed ${filterdir}/final_chimp_human_hybrid_replace_coord_get_seq.txt -fo ${filterdir}/final_chimp_human_hybrid_seq.fa

rm -f ${filterdir}/final_chimp_human_hybrid_replace_coord_get_seq.txt

######## get maf Coordinates ###########

awk '{print $7"\t"$8"\t"$9"\t"$7"|"$8"|"$9}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1 -k2,2n | uniq  > ${filterdir}/cleaned_deletion_overlapping_conserved_conserved_coord_only.txt
awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq  > ${filterdir}/cleaned_conserved_overlapping_deletion_conserved_coord_only.txt

mkdir -p ${filterdir}/maf
#get the MAFs in the coordinate files
#for conserved sequences overlapping deletions
mafsInRegion -outDir ${filterdir}/cleaned_conserved_overlapping_deletion_conserved_coord_only.txt ${filterdir}/maf ${refdir}/target_panTro4_multiz.maf

#for deleted sequences overlapping conserved elements
mafsInRegion -outDir ${filterdir}/cleaned_deletion_overlapping_conserved_conserved_coord_only.txt ${filterdir}/maf ${refdir}/target_panTro4_multiz.maf

#from the VCFs check to see how many of the fixed indels overlap with the original conserved regions, then check to see if they
#overlap with original configurations
#original configurations == original places where pairwise alignment claimed deletion was
#create script which checks that MSA sequence is present across all species given bed coordinate, bed coordinate should have name of MAF file in its 4th field
#check to see the organisms that align to it (is the sequence present in all primates? - only in a subset?)
#seperate genomes for large deletions? deletions overlapping conserved regions

#check to see if the conserved sequence of interest 
#refdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way
#mafsInRegion -outDir ${filterdir}/maf ${filterdir}/cleaned_conserved_overlapping_deletion_conserved_coord_only.txt ${refdir}/target_panTro4_multiz.maf 

######## statistics ######

#number of deleted region under each conserved region, counts in each category
#awk '{print $1"|"$2"|"$3}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort | uniq -c | awk '{print $1}' | sort -k1,1n | uniq -c 
#number of conserved region under each deletion region, counts in each category
#awk '{print $1"|"$2"|"$3}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt  | sort | uniq -c | awk '{print $1}' | sort -k1,1n | uniq -c 

#total number of deleted sites
#awk '{print $1"|"$2"|"$3}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt  | sort | uniq | wc -l

#total number of conserved sites
#awk '{print $1"|"$2"|"$3}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort | uniq | wc -l

#check to see if any human region is sustantially larger than deleted region  
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1 -k2,2n | uniq | awk '$NF>$(NF-1) {print}' | wc -l

#check to see if any human region is sustantially of different size magnitude than conserved region 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq | awk '($NF-$(NF-1)> 200) || ($(NF-1)-$NF> 200) {print}' | wc -l

#check to see if any human region is substantially larger than conserved region 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq | awk '$NF-$(NF-1)> 0 {print}' | wc -l

############ possible additional filters ? #################

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq | awk '$NF-$(NF-1)> 10 {print}' | wc -l
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1 -k2,2n | uniq | awk '$NF-$(NF-1)> 10 {print}' | wc -l

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq | awk '$NF/$(NF-1)> 1.05 && $NF/$(NF-1)<1.1 {print}' | sort -k7,7n | tail


awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq | awk '$NF/$(NF-1)> 1.05 {print}' | wc -l
#156
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$3-$2"\t"$6-$5}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt | sort -k1,1 -k2,2n | uniq | awk '$NF/$(NF-1)> 1.05 {print}' | wc -l
#42 for non


#in creating the chimp/human reference genome, I should tile sites so that each site is seperated by a minimum of x distance, this is to remove interference if other sites
#are very close to each other
#sort human regions by bed coordinate, then partition each sequence to a new reference genome if there is overlap of x distance
#extend coordinates to see how many human regions are close to each other - this should be a consequence from algorithm before
#I don't have to worry about multiple conserved sites overlapping deleted elements
#########################

join -v1 -14 -25 -o'1.1 1.2 1.3 1.4' <(sort -k4,4 ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only.bed ) <(sort -k5,5 ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised.bed) | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_conserved_only_no_partial.bed
#do the partials overlap with anything else? yes they do
#check to see if the deletion lies completely in the partial sequences 
#if not, the partials are partials themselves, and the algorithm needs to be rerun
bedtools intersect -loj -a ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed | awk '$7!=-1{print }' | awk '{print $3-$2}'

bedtools intersect -loj -a ${filterdir}/chimp_partial_deletion_coordinates_cleaned_revised.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_cleaned_deletion_only.bed |  awk '$7!=-1{print }' | awk '!(((($2>$7) && ($3<$8)) ||  (($2==$7) && ($3<$8)) || (($2>$7) && ($3==$8))) || (($2<=$7) && ($3>=$8))) {print}'

#awk '{print $0"\t"$NF-$(NF-1)}'  ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed | sort -k13,13n | tail


#do we insert deleted site before or after the human specified base?, should be after
#for partial sequences 
#probably need to liftover partial sequences, or not since the sequence cutoff is where the partial sequence ends
#end point and start point remain the same
#########################

#need a script that converts the deletion coordinates to chimp/human hybrid genome coordinates later, 
#this will be able to assess/check if the chimp annotated sites are truly polymorphic
#do I need to liftover the partial deletions back to the original human coordinates?
#the importance of partial deletions is to 


########as a check, the lift over coordinates should match######
awk '{print $1"\t"$2"\t"$3"\t"$NF}' ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed > ${filterdir}/test.txt
chimp_human_pairwise_path=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38
liftOver -minMatch=0.0001 ${filterdir}/test.txt ${chimp_human_pairwise_path}/panTro4ToHg38.over.chain.gz ${filterdir}/test_mapped.txt ${filterdir}/test_notmapped.txt
head ${filterdir}/test_mapped.txt
awk '{print $1"|"$2"|"$3"\t"$4}' ${filterdir}/test_mapped.txt | awk '$1!=$2{print}' | head
#################################################################




#for partials use this algorithm:
#the deleted conserved element belongs to the deleted one, the other flanking sequence gets all the sequence outside the deleted element
#the human mapped sequence is the same for both deletion and flanking sequence

#I have these files in the end:
#for deletions, need to keep track of conserved site - thats where the variant call matters
#what about deletions overlapping multiple conserved elements?, how to keep track of those?, each line should have one deletion, along with a meta column of all the conserved sites removed
#similar for conserved swap only - each line should have one conserved element, along with subsequent lines being deletions
#also need to filter for human gap being much larger than chimp deletion gap

#deletions which fully overlap conserved area can also partially overlap another conserved area
#so combine both deletion sites to create meta file

#similarily, conserved regions which fully overlap deletions regions may also partially overlap conserved areas
#thus, combine these files into one category as well
awk '{if(metainfo[$1"\t"$2"\t"$3"\t"$4]==""){ metainfo[$1"\t"$2"\t"$3"\t"$4]=$5"|"$6"|"$7"|"$8 } else{ metainfo[$1"\t"$2"\t"$3"\t"$4] = metainfo[$1"\t"$2"\t"$3"\t"$4]":"$5"|"$6"|"$7"|"$8 } } END{ for (id in metainfo) { print id"\t"metainfo[id]} }' ${filterdir}/chimp_deletion_overlapping_conserved_coordinates.bed | head


##conserved site completely overlapping deletions that does not intersect with larger deletion
head ${filterdir}/chimp_deletion_mostconserved_overlapping_deleted_cleaned.bed
bedtools intersect -a ${filterdir}/chimp_deletion_mostconserved_overlapping_deleted_cleaned.bed -b  

##deletion completely overlapping conserved site
head ${filterdir}/chimp_deletion_overlapping_conserved_coordinates.bed

#human gap larger than chimp gap
awk -F'[\t|]' '$10>$7{print $10-$7"\t"$0}' ${filterdir}/chimp_deletion_overlapping_conserved_coordinates.bed | sort -k1,1n | tail

##deletion partially overlapping conserved site
head ${filterdir}/chimp_partial_deletion_coordinates.bed

##conserved site paritially overlapping conserved site (revised coordinates)
head ${filterdir}/chimp_partial_deletion_coordinates_cleaned.bed

#create hybrid genome from these files

############### MISC ####################


#get the human coordinates from most conserved deleted sites
#awk -F '\t' '{print $4}' ${filterdir}/chimp_deletion_mostconserved_coordinates.bed | awk -F '|' '{print $5"\t"$6"\t"$6+$7"\t"$0}' | sort -k1,1 -k2,2n | uniq > ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only.bed


#remove deletions which remove the same human sequence (otherwise, don't know which sequence to swap in human swap) - these are repetitive sequences in chimp which removed a corresponding human sequence?
#for those deletions, its not possible to find a precise breakpoint 
#only two sequences removed 
#chr1	87604124	87604124	chr1_GL389237_random|31661|31856|195|chr1|87604124|0|+
#chr1	87604124	87604124	chr1|88387999|88388194|195|chr1|87604124|0|+
#bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only.bed -b ${filterdir}/chimp_deletion_mostconserved_coordinates_human_only.bed |  awk '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n | uniq -c  |  awk '$1>=1{print $2"\t"$3"\t"$4"\t"$5}' > ${filterdir}/chimp_deletion_mostconserved_nonoverlapping_deletions_coordinates_human_only.txt



awk -F'\t' '{print $5"\t"$6"\t"$7"\t"$8}' ${filterdir}/chimp_deletion_mostconserved_coordinates.bed | awk -F'[\t|]' '{print $4"\t"$5"\t"$6"\t"$1"|"$2"|"$3}' > ${filterdir}/chimp_deletion_mostconserved_liftedover_humand_coordinates.bed


##final list for partial deleted sites##
#this list will be the list where human sequences intersect each other 
#for these, combine the human sequence from deletion as well as human sequence from flanking region
#need information such as exactly where human breakpoints of deletion lies, as well as the extension from the flanking region
#probably need a python script to parse that 

#these are human coordinates
filterdir=/cluster_path/ape_project/deletions_project/filtering
bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_nonoverlapping_deletions_coordinates_partial_fully_human_only.txt -b ${filterdir}/chimp_deletion_mostconserved_liftedover_humand_coordinates.bed | awk '(($2==$6) || ($3==$7) || $6==-1){print}' | awk '{print $4}' | awk -F'[\t|]' '{print $1"\t"$2"\t"$3"\t"$0 }' | sort | uniq> ${filterdir}/chimp_deletion_mostconserved_nonoverlapping_deletions_coordinates_partial_fully_cleaned.txt 
filterdir=/cluster_path/ape_project/deletions_project/filtering_new

bedtools intersect -loj -a ${filterdir}/chimp_deletion_mostconserved_nonoverlapping_deletions_coordinates_partial_fully_human_only.txt -b ${filterdir}/chimp_deletion_mostconserved_liftedover_humand_coordinates.bed | awk '(($2==$6) || ($3==$7)){print}' | awk '{print $4}' 

