#this file extracts the polymorphic indels and obtains info on the allele frequency, the populations that they are varying, etc.
#it should be run after analyze_vcf_stats_submit.sh
genome_file_id=$1
af_filter=$2
sgdp_output_dir=$3   
sgdp_meta_file=$4
vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files

del_dir=/cluster_path/ape_project/deletions_project 

 
rm -f ${vcf_files}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt

#B stands for biallelic, NB stands for not biallelic, in the end we can filter to keep only biallelic variants
vcftools --gzvcf ${vcf_files}/${genome_file_id}.flt.vcf.gz --keep-only-indels --non-ref-af ${af_filter} --counts --out ${vcf_files}/${genome_file_id}_af_${af_filter}
awk -F"[\t:]" '(NR>1) && ($4!=0){if( NF> 8) {print $1"_"$2-1"\t"$6"\t"$4"\t"$6/$4"\tNB"} else{print $1"_"$2-1"\t"$6"\t"$4"\t"$6/$4"\tB" } }' ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.count > ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info


grep "#CHROM" ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf | sed "s#/cluster_path_temp/sgdp/${genome_file_id}_bam_files/##g" | sed "s#/broad/hptmp/sreilly/James/${genome_file_id}_bam_files/##g" | sed "s#_${genome_file_id}.srt.bam##g" | awk '{for(i=10;i<=NF;i++){print i-9"\t"$i}}' | sort -k2,2 > ${vcf_files}/${genome_file_id}_af_${af_filter}_column_number_pop_map.txt


join -12 ${vcf_files}/${genome_file_id}_af_${af_filter}_column_number_pop_map.txt  -21 ${sgdp_meta_file} | awk '{print $2"\t"$1"\t"$3"\t"$4}' | sort -k1,1 > ${vcf_files}/${genome_file_id}_column_number_pop_map_full_info.txt 

#variant ID will be chimp del seq, human coord del seq, inv info - this ID will be used to merge with the files from check_deletion_validity_maf
#file inputted into script will be variantID, VCFmatchID, then each of the allele info entries from the VCF
inputFile=${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt
numCol=$(awk '{print NF}' ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt | head -1) 
paste <(cut -f13-19 ${inputFile} |  awk 'BEGIN{OFS="|";} {$1=$1}1' | awk '{print "C""\t"$0}' )  <(cut -f20-21 ${inputFile} |  awk 'BEGIN{OFS="_";} {$1=$1}1' ) <( awk '{ for(i=30;i<=NF-2;i++){printf $i"\t"};printf $(NF-1)"\n" }' ${inputFile} ) > ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_parsed.txt

inputFile=${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt 
numCol=$(awk '{print NF}' ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt  | head -1)
paste <(cut -f13-19 ${inputFile} |  awk 'BEGIN{OFS="|";} {$1=$1}1' | awk '{print "P""\t"$0}'  ) <(cut -f20-21 ${inputFile} | awk 'BEGIN{OFS="_";} {$1=$1}1') <( awk '{ for(i=30;i<=NF-2;i++){printf $i"\t"};printf $(NF-1)"\n" }' ${inputFile}  ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_parsed.txt


#remember that the first 6 coordinates of pass align are the chimp/human hybrid coordinates
#the next 3 afterwards are the coordinates that completely overlie conserved seq/chimp seq.
#so for deletion overlapping conserved it is the deletion seq
inputFile=${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt 
numCol=$(awk '{print NF}' ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | head -1)
paste <(cut -f7-12,19 ${inputFile} |  awk 'BEGIN{OFS="|";} {$1=$1}1' | awk '{print "D""\t"$0}'  )  <(cut -f20-21 ${inputFile} | awk 'BEGIN{OFS="_";} {$1=$1}1' ) <( awk '{ for(i=30;i<=NF-2;i++){printf $i"\t"};printf $(NF-1)"\n" }' ${inputFile}  ) >> ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_parsed.txt


#now run python script to calculate the population info
python extractPopInfoFromVCF.py ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_parsed.txt ${vcf_files}/${genome_file_id}_column_number_pop_map_full_info.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt 4

#rearrange output file
paste <(awk '{print $3"\t"$2"\t"$1}' ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt ) <(cut -f4- ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt ) > ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info_revised.txt

mv ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info_revised.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt 

#sort to use join later
sort -k1,1 ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt  > ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info_sorted.txt 
mv ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info_sorted.txt ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt

sort -k1,1 ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info > ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info.sorted
mv ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info.sorted ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info

#combine with actual allele frequency info 
join -11 ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_cleaned_info.txt -21 ${vcf_files}/${genome_file_id}_af_${af_filter}.frq.info > ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_with_freq.txt




