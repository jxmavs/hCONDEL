genome_file_id=$1
af_filter=$2
sgdp_output_dir=$3

#genome_file_id=ChimpHumanGenome_Partition_2
#af_filter=1
#sgdp_output_dir=/cluster_path_temp/sgdp

#This script performs the initial filter
#the initial filter used is to only keep deletion sites which are deleted in all humans and are correctly annotated to remove the same sequence

#This script produces general statistics on the sites after aligning
#but is also important to filter out sites you want from the ones you don't


#this script looks at the vcf files and analyzes basic statistics
temp_dir=/cluster_path_temp/sgdp/temp
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
#input_vcf_file=${genome_file_id}.flt.vcf.gz
rm -f ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

#change af_filter to relax complete deletion definition
vcftools --gzvcf ${vcf_files}/${genome_file_id}.flt.vcf.gz --keep-only-indels --non-ref-af ${af_filter} --recode --out ${vcf_files}/${genome_file_id}_af_${af_filter}
vcf2bed --deletions < ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf | awk '$8=="PASS" {print}' >  ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed

#$(NF-5)"\t"$(NF-4)"\t"$(NF-3) are the converted deletion coordinates
#$(NF-2)"\t"$(NF-1)"\t"$NF are the converted conserved coordinates
count=$( awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of conserved sites (conserved overlapping deletion) total" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $(NF-2)"\t"$(NF-1)"\t"$NF}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of deleted sites (conserved overlapping deletion) total" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt
echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of deleted sites (deletion overlapping conserved) total" ) <( echo ${count} )  >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $(NF-2)"\t"$(NF-1)"\t"$NF}'  ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of conserved sites in deleted sites (deletion overlapping conserved) total" ) <( echo ${count} )  >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt
echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of partial conserved sites (partial conserved overlapping deletion) total" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $(NF-2)"\t"$(NF-1)"\t"$NF}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort | uniq | wc -l | awk '{print $1}' )
paste <( echo "prelim number of deleted sites (partial conserved overlapping deletion) total" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt
echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

################# how many lie completely in conserved region (conserved sequences completely overlapping deletions)  #################################
#$(NF-5)"\t"$(NF-4)"\t"$(NF-3) is the new converted conserved coordinate

awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt

num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt | head -1 ) 

#intersect conserved coordinates overlapping deletion sites with vcf file of allele frequency of 1 (complete deletions)
bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)!=-1 {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion.txt

#need to add 1 to coordinate since vcf takes in bp right before deletion
count=$( awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print $1"\t"$2"\t"$3 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion.txt | sort -k1,1 -k2,2n | uniq  | wc -l | awk '{print $1}' ) 
paste <( echo "number conserved overlapping deletion" ) <( echo ${count} )>> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of deletion size from conserved overlapping deletion:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print $4"\t"$(num_fields+2)+1"\t"$(num_fields+3) }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion.txt | sort -k1,1 -k2,2n | uniq   | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt


######### how many remove exact sequence of conserved overlapping deletion? #########################################
#new coordinates, new coordinates of gap,  original conserved coordinates, human conserved coord, original gap coord, human gap coord

awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_revised_${genome_file_id}.txt

extendBP=1000

#extend coordinates to get more context in the check
python extendCoordinatesV2.py ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_revised_${genome_file_id}.txt ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.sizes ${extendBP} ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt 


#redo intersect to pick up deletions possibly outside of region, we do this because of the left-alignment issue
num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt | head -1 )

#only look at positions where the entire VCF annotated deletion lies within the extended coordinates
bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt  -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_extended.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_extended.txt | sort | uniq > ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_conserved_coord_extended.txt

bedtools getfasta -name -fi ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa -bed ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_conserved_coord_extended.txt -fo ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_conserved_coord_extended.fa

#the output is the same bed intersected file with an extra field of 1 if it passed the filtering, or 0 if not
python deletionCheckOverlapV3.py  ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_extended.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_conserved_coord_extended.fa 1 7 $((${num_fields}+1)) ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_extended_checked.txt

#for now only take sequences that remove same sequence, ignore the extended coordinates
#cut -f4- is to remove the extended coordinates
#THIS IS ALSO WHERE WE KEEP THE FILES FOR FURTHER MAF FILTERING IN CHECK_DELETION_VALIDITY.SH script
awk '$NF==1{print}' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_extended_checked.txt | cut -f4- > ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt

#################

#7,8,9 is the original conserved coordinates
count=$(awk '{print $7"\t"$8"\t"$9 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt  | sort -k1,1 -k2,2n | uniq | wc -l)
paste <(echo "number of unique conserved regions - conserved overlapping deletion with gap coordinates matching") <(echo ${count}) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

#13,14,15 is the original deletion coordinates lying under conserved
count=$(awk '{print $13"\t"$14"\t"$15 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt  | sort -k1,1 -k2,2n | uniq | wc -l)
paste <(echo "number of unique deletion regions - conserved overlapping deletion with gap coordinates matching") <(echo ${count}) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt


echo -e "distribution of conserved region sizes - conserved overlapping deletion with deleted sequences matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk  '{print $7"\t"$8"\t"$9 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt | awk '{print}' | sort | uniq | sort -k1,1 -k2,2n | uniq  | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of deletion region sizes - conserved overlapping deletion with deleted sequences matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk  '{print $13"\t"$14"\t"$15 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt | awk '{print}' | sort | uniq  | sort -k1,1 -k2,2n | uniq  | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

#lie partial - good to know
#awk '{print $1"\t"$2"\t"$3 }' ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.txt | sort -k1,1 -k2,2n | uniq  | wc -l

#bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt -b ${vcf_files}/${input_vcf_file}_af_${af_filter}.recode.vcf.bed -loj  | awk '$6!=-1 {print $1"\t"$2"\t"$3}' | sort | uniq | wc -l
#bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj  | awk '$6!=-1 {print $1"\t"$2"\t"$3}' | sort | uniq | wc -l


################# get the conserved coordinates lying underneath the deletion site ##############################

#first look at the small sites, then look at the large sites (structural variant deletions)
#in this section, only look at the initial number of deletions annotated as lying over conserved sites
#cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt is different than ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt

num_fields=$( awk '{print NF}'  ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt  | head -1 )

#new coordinates of conserved element  original deletion coordinates, human deletion coord, original gap coord, human gap coord
awk -v num_fields=${num_fields} '{print $(NF-2)"\t"$(NF-1)"\t"$NF"\t"$0}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt | cut -f1-${num_fields} |  sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_under_only_${genome_file_id}.txt

num_fields=$( awk '{print NF}'  ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_under_only_${genome_file_id}.txt | head -1 )

bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_under_only_${genome_file_id}.txt -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)!=-1' > ${vcf_files}/${genome_file_id}_af_${af_filter}.underlie.deletion.conserved.txt

# deletion overarching count #

count=$( awk -v num_fields=${num_fields} '{print $(num_fields+1)"\t"$(num_fields+2)"\t"$(num_fields+3) }'  ${vcf_files}/${genome_file_id}_af_${af_filter}.underlie.deletion.conserved.txt | sort | uniq | wc -l  | awk '{print $1}' )
paste <( echo "number of deletions overlapping conserved" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

############### deletion overlappiing conserved exact coordinate overlap #######################

#new coordinates of deletion, new coordinates of conserved element,  original deletion coordinates, human deletion coord, original conserved coord, human conserved coord
awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_revised_${genome_file_id}.txt

extendBP=1000

#extend coordinates to get more context in the check
python extendCoordinatesV2.py ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_revised_${genome_file_id}.txt ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.sizes ${extendBP} ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_revised_${genome_file_id}_extended.txt 

#redo intersect to pick up conserveds possibly outside of region, we do this because of the left-alignment issue
num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_revised_${genome_file_id}_extended.txt | head -1 )

bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_revised_${genome_file_id}_extended.txt  -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3  {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_extended.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_extended.txt | sort | uniq > ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_deletion_coord_extended.txt

bedtools getfasta -name -fi ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa -bed ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_deletion_coord_extended.txt -fo ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_deletion_coord_extended.fa

#the output is the same bed intersected file with an extra field of 1 if it passed the filtering, or 0 if not
python deletionCheckOverlapV3.py  ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_extended.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlapping_conserved_deletion_coord_extended.fa 1 7 $((${num_fields}+1)) ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_extended_checked.txt

#for now only take sequences that remove same sequence, ignore the extended coordinates
#THIS IS ALSO WHERE WE KEEP THE FILES FOR FURTHER MAF FILTERING IN CHECK_DELETION_VALIDITY.SH script
awk '$NF==1{print}' ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_extended_checked.txt | cut -f4- > ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt


#count=$( awk -v num_fields=${num_fields} '$(num_fields-1)==$(num_fields+2)+1 && $(num_fields)==$(num_fields+3){print }' ${vcf_files}/${genome_file_id}_af_${af_filter}.underlie.deletion.conserved.txt | awk -v num_fields=${num_fields} '{print $(num_fields-2)"\t"$(num_fields-1)"\t"$(num_fields)}' | sort | uniq | wc -l )

count=$( awk '{print $7"\t"$8"\t"$9 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | sort | uniq | wc -l )

paste <( echo "number of unique deletion regions - deletions overlapping conserved exact coordinate match" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

count=$( awk '{print $13"\t"$14"\t"$15 }'  ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | sort | uniq | wc -l )

paste <( echo "number of unique conserved regions - deletions overlapping conserved exact coordinate match" ) <( echo ${count} ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of deletion region sizes - deletion overlapping conserved with with deleted sequencees matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk '{print $7"\t"$8"\t"$9 }'  ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | awk '{print}' | sort | uniq | awk '{print $3-$2}' | sort -k1,1n | uniq -c >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of conserved region sizes - deletion overlapping conserved with with deleted sequences matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk '{print $13"\t"$14"\t"$15 }'  ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | awk '{print}' | sort | uniq| awk '{print $3-$2}' | sort -k1,1n | uniq -c >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

############ partial check - pretty much the same as conserved over deletion #####################################################


awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt

num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt | head -1 ) 

#intersect partial conserved coordinates overlapping deletion sites with vcf file of allele frequency of 1 (complete deletions)
bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_over_only_${genome_file_id}.txt -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)!=-1 {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion.txt

#need to add 1 to coordinate since vcf takes in bp right before deletion
count=$( awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print $1"\t"$2"\t"$3 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion.txt | sort -k1,1 -k2,2n | uniq  | wc -l | awk '{print $1}' ) 
paste <( echo "number partial conserved overlapping deletion" ) <( echo ${count} )>> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of deletion size from partial_conserved overlapping deletion:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print $4"\t"$(num_fields+2)+1"\t"$(num_fields+3) }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion.txt | sort -k1,1 -k2,2n | uniq   | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt


######### how many remove exact sequence of partial_conserved overlapping deletion? #########################################
#new coordinates, new coordinates of gap,  original partial_conserved coordinates, human partial_conserved coord, original gap coord, human gap coord

awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_revised_${genome_file_id}.txt

extendBP=1000

#extend coordinates to get more context in the check
python extendCoordinatesV2.py ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_revised_${genome_file_id}.txt ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.sizes ${extendBP} ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt 


#redo intersect to pick up deletions possibly outside of region, we do this because of the left-alignment issue
num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt | head -1 )

#only look at positions where the entire VCF annotated deletion lies within the extended coordinates
bedtools intersect -sorted -a ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_revised_${genome_file_id}_extended.txt  -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_extended.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_extended.txt | sort | uniq > ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_partial_conserved_coord_extended.txt

bedtools getfasta -name -fi ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa -bed ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_partial_conserved_coord_extended.txt -fo ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_partial_conserved_coord_extended.fa

#the output is the same bed intersected file with an extra field of 1 if it passed the filtering, or 0 if not
python deletionCheckOverlapV3.py  ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_extended.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlapping_deletion_partial_conserved_coord_extended.fa 1 7 $((${num_fields}+1)) ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_extended_checked.txt

#for now only take sequences that remove same sequence, ignore the extended coordinates
#THIS IS ALSO WHERE WE KEEP THE FILES FOR FURTHER MAF FILTERING IN CHECK_DELETION_VALIDITY.SH script
awk '$NF==1{print}' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_extended_checked.txt | cut -f4- > ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt

#################

#7,8,9 is the original partial_conserved coordinates
count=$(awk '{print $7"\t"$8"\t"$9 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt  | sort -k1,1 -k2,2n | uniq | wc -l)
paste <(echo "number of unique partial_conserved regions - partial_conserved overlapping deletion with gap coordinates matching") <(echo ${count}) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

#13,14,15 is the original deletion coordinates lying under partial_conserved
count=$(awk '{print $13"\t"$14"\t"$15 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt  | sort -k1,1 -k2,2n | uniq | wc -l)
paste <(echo "number of unique deletion regions - partial_conserved overlapping deletion with gap coordinates matching") <(echo ${count}) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt


echo -e "distribution of partial conserved region sizes - partial_conserved overlapping deletion with deleted sequences matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk  '{print $7"\t"$8"\t"$9 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt | awk '{print}' | sort | uniq | sort -k1,1 -k2,2n | uniq  | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "distribution of deletion region sizes - partial_conserved overlapping deletion with deleted sequences matching:\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

awk  '{print $13"\t"$14"\t"$15 }' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt | awk '{print}' | sort | uniq  | sort -k1,1 -k2,2n | uniq  | awk '{print $3-$2}' | sort -k1,1n | uniq -c >>  ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

echo -e "\n" >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt


##### combine all coordinates together to check for deletion validity #############
#${vcf_files}/${genome_file_id}_af_${af_filter}.underlie.deletion.conserved.txt

#structural variants look through each individual person file
#ls -l /broad/hptmp/sreilly/James/${genome_file_id}_vcf_files/*.sv.vcf.gz | awk -F'_' | '{print $1"_"$2"\t"$0}' >  ${temp_dir}/sv_file_list_${genome_file_id}.txt
ls -l ${sgdp_output_dir}/${genome_file_id}_vcf_files/*.sv.vcf.gz | awk '{print $NF}' | awk -F'/' '{print $NF}' | awk -F'_ChimpHumanGenome_' '{print $1"\t"$0}' > ${temp_dir}/sv_file_list_${genome_file_id}.txt

num_samples=$( wc -l ${temp_dir}/sv_file_list_${genome_file_id}.txt | awk '{print $1}' )

num_fields=$( awk '{print NF}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_under_only_${genome_file_id}.txt | head -1)

rm -f  ${vcf_files}/all_sv_${genome_file_id}_overlapped_sites.txt
while read -a array
do
    sample_id=${array[0]}
    sv_file=${array[1]}
    gunzip -c ${vcf_files}/${sv_file} > ${vcf_files}/${sample_id}.sv.vcf
    
    bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_under_only_${genome_file_id}.txt -b ${vcf_files}/${sample_id}.sv.vcf -loj | awk -v num_fields=${num_fields} '$(num_fields+2)!=-1 && $(num_fields+5)=="<DEL>" {print}' > ${vcf_files}/${sample_id}_af_${af_filter}.underlie.deletion.conserved_sv_intersected.txt
    #use python script to assess exact deletion?
 
    #remove the first 3 coordinates which are the chimp/human hybrid genome coord used for intersection
    awk  -v sample_id=${sample_id} '{print $0"\t"sample_id}' ${vcf_files}/${sample_id}_af_${af_filter}.underlie.deletion.conserved_sv_intersected.txt | cut -f4-${num_fields}  | sort | uniq >> ${vcf_files}/all_sv_${genome_file_id}_overlapped_sites.txt
    
    rm -f ${vcf_files}/${sample_id}.sv.vcf

done<${temp_dir}/sv_file_list_${genome_file_id}.txt

#count the number of sites deleted in all humans
#THIS IS ALSO WHERE WE KEEP THE FILES FOR FURTHER MAF FILTERING IN CHECK_DELETION_VALIDITY.SH script
awk '{print $0}' ${vcf_files}/all_sv_${genome_file_id}_overlapped_sites.txt | sort | uniq -c | awk -v num_samples=${num_samples} '$1==num_samples{print}'| cut -f2- > ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_SV_pass_align.txt


count=$( wc -l ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_SV_pass_align.txt | awk '{print $1}' )
paste <(echo "number of structural variants deleted in all humans:") <(echo ${count}) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_vcf_stats.txt

#keep the elements where deletion lies completely over conserved element

