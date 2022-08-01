#continuation of before after creation of chimp/human hybrid genome
#after running finalFilteringProtocol_v3.sh, run the following below to create the chimp/human hybrid genome
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp
#seperate out coordinates here:
prefixFileName="ChimpHumanGenome_Partition"
####### be careful with the following code below!!!, don't remove files that you need #########
rm -f /cluster_path/ape_project/deletions_project/chimp_human_hybrid/${prefixFileName}* 
bash seperateOutChimpHumanHybridCoordinates.sh ${prefixFileName}

######### convert each of the genome coordinates to fasta ############## 
chimp_human_ref_genome_file_names_file=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
ls -l /cluster_path/ape_project/deletions_project/chimp_human_hybrid/${prefixFileName}*.txt | awk '{print $NF}' | awk -F'[/.]' '{print $(NF-1)}' > ${chimp_human_ref_genome_file_names_file}

chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid
#hybrid_run_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt

awk '{print $1"|"$2"|"$3"\t"$0}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_with_id.txt
awk '{print $1"|"$2"|"$3"\t"$0}' ${filterdir}/cleaned_deletion_overlapping_conserved.txt  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_with_id.txt
awk '{print $1"|"$2"|"$3"\t"$0}' ${filterdir}/cleaned_partial_conserved_overlapping_deletion.txt  > ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_with_id.txt

#rm -f ${hybrid_run_file_names}
#preconverted file is simply in the same format as the original version, except with only the relevant columns pulled out for constructing the hybrid genome
while read -a array
do
	file_name=${array[0]}
	#get the fasta file sequences for the particular genome iteration
	awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3"|"$NF}' ${chimp_human_hybrid_dir}/${file_name}.txt > ${chimp_human_hybrid_dir}/${file_name}_get_seq_temp.txt
	bedtools getfasta -name -fi ${chimp_ref_dir}/panTro4.fa -bed ${chimp_human_hybrid_dir}/${file_name}_get_seq_temp.txt -fo ${chimp_human_hybrid_dir}/${file_name}.fa
	
	rm -f ${chimp_human_hybrid_dir}/${file_name}_get_seq_temp.txt
	#extract the other coordinates into the seperate files
	awk '{print $1"|"$2"|"$3}' ${chimp_human_hybrid_dir}/${file_name}.txt > ${chimp_human_hybrid_dir}/${file_name}_with_ids.txt
	
	#get the conserved overlapping deleted elements for the particular hybrid genome
	join -11 -21 <(sort -k1,1 ${chimp_human_hybrid_dir}/${file_name}_with_ids.txt) <(sort ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_with_id.txt )| awk '{for(i=2;i<=NF;i++)printf $i"\t"; printf"\n"}' >  ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_preconverted_${file_name}.txt

	#get the deleted overlapping conserved elements for the particular hybrid genome
	join -11 <(sort -k1,1 ${chimp_human_hybrid_dir}/${file_name}_with_ids.txt) -21 <(sort ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_with_id.txt )| awk '{for(i=2;i<=NF;i++)printf $i"\t"; printf"\n"}' >  ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_preconverted_${file_name}.txt
	
	#get the deleted overlapping conserved elements for the particular hybrid genome
	join -11 <(sort -k1,1 ${chimp_human_hybrid_dir}/${file_name}_with_ids.txt) -21 <(sort ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_with_id.txt)| awk '{for(i=2;i<=NF;i++)printf $i"\t"; printf"\n"}' >  ${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_preconverted_${file_name}.txt
	
	#echo -e "${file_name}\t${chimp_human_hybrid_dir}/${file_name}.txt\t${chimp_human_hybrid_dir}/${file_name}.fa\t${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.txt\t${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.fa" >> ${hybrid_run_file_names}
	
	#${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.txt
	#${chimp_human_hybrid_dir}/${file_name}_hybrid_out_coord.fa
	#above are the new coordinates and fasta file after creating the chimp/human hybrid genome

done<${chimp_human_ref_genome_file_names_file}
############################################################

#these other scripts mainly take info from this file:
#/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
#:%s/all_runs/all_runs_convert_to_hybrid/g
#:%s/submit_all/submit_convert_to_hybrid/g
bash create_hybrid_sequence_submit.sh

#get the sizes of the hybrid sequence after creating them
bash get_hybrid_sequence_sizes_submit.sh

#get bwa indeces 
bash build_chimp_hybrid_index_submit.sh


#then run the following to convert the underlying coordinates in each region to the coordinates corresponding to the hybrid genome
bash convert_to_chimp_human_hybrid_submit.sh 
#converted coordinates are at the end, the last 6 coordinates
#from the last 6 fields, the first 3 are the overlapping coordinates
#the last 3 are the coordinate lying underneath the overlapping coordinates

#### below are the computationally intensive parts ####

#now align the unitigs from simons genome diversity project to chimp/human hybrid
#below calls script bwa_align_tigs.sh
#cd /cluster_path_temp/sgdp
#ls -l | awk 'NR>1 {print $9}' | awk -F'.'  '{print $0"\t"$1}' | tail -10 | head -5 >sgdp_file_names_3.txt
#ls -l *.mag.gz | awk '{print $9}' | awk -F'.'  '{print $0"\t"$1}' | sort >all_sgdp_file_names.txt

#compare with dataset advertised
#wget http://simonsfoundation.s3.amazonaws.com/share/SCDA/datasets/10_24_2014_SGDP_metainformation_update.txt
#sed 's/\r$//' 10_24_2014_SGDP_metainformation_update.txt 
#awk '(NR>1) && ($1=="C"){print $2}' 10_24_2014_SGDP_metainformation_update.txt| sort > sgdp_real_dataset.txt
#awk '{print $2}' all_sgdp_file_names.txt > all_sgdp_file_names_names_only.txt
#diff all_sgdp_file_names_names_only.txt sgdp_real_dataset.txt
#awk 'NR==FNR{a[$1]=1;next} ($1 in a) {print}'  sgdp_data_paper_all.txt all_sgdp_file_names_names_only.txt 
#get missing ids
#awk 'NR==FNR{a[$1]=1;next} !($1 in a) {print}'  all_sgdp_file_names_names_only.txt sgdp_paper_all_with_embargo_status.txt | awk '$3=="X"{print $1}' | sort > sgdp_missing_ids.txt
#check embargo status - make sure that they are all X - open access
#awk 'NR==FNR{a[$1]=1;next} ($1 in a) {print}'  all_sgdp_file_names_names_only.txt sgdp_paper_all_with_embargo_status.txt | awk '{print $3}' | sort | uniq
#del

awk 'NR<=2{print $1}' /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt > /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1_2.txt 

awk 'NR>=3{print $1}' /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt > /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3_4.txt 
chmod 755 /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1_2.txt 
chmod 755 /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3_4.txt 

sgdp_input_dir=/cluster_path_temp/sgdp
#sgdp_output_dir=/cluster_path_temp/sgdp
#sgdp_file_name=sgdp_file_names.txt
#sgdp_output_dir=/broad/hptmp/sreilly/James
sgdp_output_dir=/cluster_path_temp/sgdp
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3_4.txt 
#mv ChimpHumanGenome_Partition_1_bam_files /cluster_path_temp/sgdp
sgdp_file_name=all_sgdp_file_names.txt
bash align_sgdp_samples_submit.sh ${sgdp_input_dir} ${sgdp_output_dir} ${sgdp_file_name} ${chimp_ref_file_names}


#command for steve?
#bash /cluster_path/ape_project/deletions_project/align_sgdp_samples_submit.sh /cluster_path_temp/sgdp /broad/hptmp/sreilly/James all_sgdp_file_names.txt /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1_2.txt 

#check the length 
qsub get_sgdp_max_read_length.sh

#afterwards call variants, the following script sgdp_call_chimp_hybrid_variants_submit.sh runs these scripts:
#call_sgdp_SV_variants.sh, call_sgdp_small_variants.sh

awk 'NR==3{print}' /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt > /cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3.txt
sgdp_output_dir=/broad/hptmp/sreilly/James
#sgdp_output_dir=/cluster_path_temp/sgdp
sgdp_file_name=all_sgdp_file_names.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_3.txt
#chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_file_4.txt
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files_1_2.txt
bash sgdp_call_chimp_hybrid_variants_submit.sh ${sgdp_file_name} ${sgdp_output_dir} ${chimp_ref_file_names}

chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid
sgdp_output_dir=/broad/hptmp/sreilly/James
genome_file_id=ChimpHumanGenome_Partition_0
sgdp_bam_dir=${sgdp_output_dir}/${genome_file_id}_bam_files
sgdp_vcf_dir=${sgdp_output_dir}/${genome_file_id}_vcf_files
chimp_human_hybrid_out_fasta_name=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.fa
qsub call_sgdp_small_variants.sh ${chimp_human_hybrid_out_fasta_name} ${sgdp_bam_dir} ${sgdp_vcf_dir} ${genome_file_id}

#does the first pass filtering of only extracting deletions which exactly remove same site, also gets the prelim stats on all the deletions
#bash analyze_vcf_stats_submit.sh
#analyze_vcf_stats_post_call_new.sh 

#does the next pass filter on the MAF alignment, the output from that can be used to make oligos
#bash check_deletion_validity_maf.sh 

#then design the oligos - this will design the small deletion oligos, the small deletion combinatorial oligos, as well as the large deletion oligos, including designing the combinatorial ones:
#design_oligos.sh


#get the final distribution of deletion sizes and graph in R
#delGapStats.R
########## if looking at polymorphic sites ###########
#set af_filter to 0.00001 in all the scripts
bash analyze_vcf_stats_submit.sh

bash extract_freq_pop_info_submit.sh

#combine all files from extract_freq_pop_info_submit.sh
af_filter=0.00001
sgdp_output_dir=/cluster_path_temp/sgdp
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt

rm -f ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt 
while read -a array 
do
	genome_file_id=${array[0]}
	#echo ${genome_file_id}
	vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
	cat ${vcf_files}/${genome_file_id}_af_${af_filter}_all_pass_align_VCF_with_freq.txt >> ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt 
done<${chimp_ref_file_names}





bash check_deletion_validity_maf.sh 



#TODO: combine results from frequency and conservation and output into one file
paste <(cut -f7-13 ${post_vcf_filter_dir}/all_af_${af_filter}_maf_check_coord_output.txt |  awk 'BEGIN{OFS="|";} {$1=$1}1' ) <( awk '{print}' ${post_vcf_filter_dir}/all_af_${af_filter}_maf_check_coord_output.txt ) > ${post_vcf_filter_dir}/all_af_${af_filter}_maf_check_coord_output_temp.txt
join -11 ${vcf_files}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt -21 ${post_vcf_filter_dir}/all_af_${af_filter}_maf_check_coord_output_temp.txt > ${vcf_files}/all_af_${af_filter}_master_freq_cons_file.txt 

#keep only the biallelic variants
#now use R to graph 
#doubleton alt - only two chromosomes have alt
awk '$NF=="B" && $3=="C" && $(NF-2)-$(NF-3)==2{print}' ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt | awk '{print $5}' | awk -F"|" '{print $1}' | sort | uniq -c | sort -k1,1n

#doubleton ref - only two chromosomes have ref
awk '$NF=="B" && $3=="C" && $(NF-3)==2{print }' ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt | awk '{print $4}' | awk -F"|" '{print $1}' | sort | uniq -c | sort -k1,1n

#look at all variants, calculate proportion 
awk '$3=="C" {print}' ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt | wc -l

filterdir=/cluster_path/ape_project/deletions_project/filtering_new

wc -l ${filterdir}/cleaned_conserved_overlapping_deletion.txt 

################## for aligning sgdp to hg38 #######################

sgdp_input_dir=/cluster_path_temp/sgdp
sgdp_output_dir=/cluster_path_temp/sgdp/hg38
genome_file_id=simons_hg38
sgdp_file_name=/cluster_path_temp/sgdp/all_sgdp_file_names.txt
ref_file_path=/cluster_path/ape_project/deletions_project/reference_genome_index/human/bwa_hg38
bash align_sgdp_samples_hg38_submit.sh ${sgdp_input_dir} ${sgdp_output_dir} ${sgdp_file_name} ${ref_file_path} ${genome_file_id}

#call variants afterwards
sgdp_bam_dir=${sgdp_output_dir}/${genome_file_id}_bam_files
sgdp_vcf_dir=${sgdp_output_dir}/${genome_file_id}_vcf_files
mkdir -p ${sgdp_vcf_dir}
ref_fasta_file=/cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.fa
bash sgdp_call_variants_hg38_submit.sh ${sgdp_bam_dir} ${sgdp_vcf_dir} ${ref_fasta_file} ${sgdp_file_name} ${genome_file_id}
echo "bash sgdp_call_variants_hg38_submit.sh ${sgdp_bam_dir} ${sgdp_vcf_dir} ${ref_fasta_file} ${sgdp_file_name} ${genome_file_id}"
#compare hg38 sgdp with the polymorphic deletion file 
sgdp_bam_dir=
bash compare_poly_hg38_aligned_hybrid_genome_aligned.sh

af_filter=0.00001
genome_file_id=simons_hg38
vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
vcftools --gzvcf ${vcf_files}/${genome_file_id}.flt.vcf.gz --keep-only-indels --non-ref-af ${af_filter} --counts --out ${vcf_files}/${genome_file_id}_af_${af_filter}

#parse out VCF ID 
#compare 1000 bp sequence in region centered on the hCONDEL site with the VCF site to know for sure
awk '{print}' ${sgdp_output_dir}/all_af_${af_filter}_all_pass_align_VCF_with_freq.txt 

#get all the VCF files that have been merged to the initial hCONDELs
rm -f ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align.txt
while read -a array
do
	vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files
	cut -f13- ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlap_deletion_pass_align.txt | awk '{print C"\t"$0}' > ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align.txt
	cut -f13- ${vcf_files}/${genome_file_id}_af_${af_filter}_partial_conserved_overlap_deletion_pass_align.txt | awk '{print P"\t"$0}' >> ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align.txt
	cut -f7-12,19- ${vcf_files}/${genome_file_id}_af_${af_filter}_deletion_overlap_conserved_pass_align.txt | awk '{print D"\t"$0}' >>  ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align.txt
done

###### extend bp, then intersect,  then check overlap #######
extendBP=1000

#extend coordinates to get more context in the check
python extendCoordinatesV2.py ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align.txt ${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.sizes ${extendBP} ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.txt

#redo intersect to pick up deletions possibly outside of region, we do this because of the left-alignment issue
num_fields=$( awk '{print NF}' ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.txt | head -1 )

#only look at positions where the entire VCF annotated deletion lies within the extended coordinates
bedtools intersect -sorted -a ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.txt  -b ${vcf_files}/${genome_file_id}_af_${af_filter}.recode.vcf.bed -loj | awk -v num_fields=${num_fields} '$(num_fields+2)>$2 && $(num_fields+3)<$3 {print}' > ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_extended.txt

awk '{print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}' ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.txt | sort | uniq > ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended_coord_only.txt

bedtools getfasta -name -fi ${ref_fasta_file} -bed ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended_coord_only.txt -fo ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.fa

#the output is the same bed intersected file with an extra field of 1 if it passed the filtering, or 0 if not
python deletionCheckOverlapV3.py  ${vcf_files}/${genome_file_id}_af_${af_filter}_conserved_overlapping_deletion_extended.txt ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended.fa 1 7 $((${num_fields}+1)) ${sgdp_output_dir}/all_vcf_af_${af_filter}_conserved_overlap_deletion_pass_align_extended_checked.txt



##########
genome_file_id=ChimpHumanGenome_Partition_0
af_filter=0.00001
sgdp_output_dir=/cluster_path_temp/sgdp
del_dir=/cluster_path/ape_project/deletions_project
sgdp_meta_file=${del_dir}/sgdp_region_country_id_cleaned.txt
vcf_files=${sgdp_output_dir}/${genome_file_id}_vcf_files

############ test #################
bash call_sgdp_SV_variants.sh /cluster_path_temp/sgdp/ChimpHumanGenome_Partition_2_bam_files /cluster_path_temp/sgdp/ChimpHumanGenome_Partition_2_vcf_files /cluster_path/ape_project/deletions_project/chimp_human_hybrid/ChimpHumanGenome_Partition_2_hybrid_out_coord.fa S_Altaian-1 ChimpHumanGenome_Partition_2
fermi_program_path=/cluster_path/bin/fermikit/fermi.kit
#call SV
sgdp_bam_dir=/broad/hptmp/sreilly/James/ChimpHumanGenome_Partition_0_bam_files 
sgdp_vcf_dir=/broad/hptmp/sreilly/James/ChimpHumanGenome_Partition_0_vcf_files

ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/ChimpHumanGenome_Partition_0_hybrid_out_coord.fa
${fermi_program_path}/htsbox abreak -cuf ${ref_file_path} ${sgdp_bam_dir}/S_Hezhen-1_ChimpHumanGenome_Partition_0.srt.bam | gzip -1 > ${sgdp_vcf_dir}/S_Hezhen-1_ChimpHumanGenome_Partition_0.sv.vcf.gz

ref_file_path=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/ChimpHumanGenome_Partition_0_hybrid_out_coord.fa

${fermi_program_path}/htsbox abreak -cuf ${ref_file_path} /cluster_path_temp/sgdp/S_Zapotec-1.srt.bam | gzip -1 > /cluster_path_temp/sgdp/S_Zapotec-1.sv.vcf.gz


########### new test ##################
vcf_files=/cluster_path_temp/sgdp/vcf_files
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

#more stringent test to see if there exists overlap
awk '{print $10"\t"$11"\t"$12}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_hybrid_converted.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_only.txt
awk '{print $10"\t"$11"\t"$12}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_hybrid_converted.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_only.txt
#awk '{print $7"\t"$8"\t"$9}' ${filterdir}/cleaned_conserved_overlapping_deletion.txt | sort -k1,1 -k2,2n | uniq  > ${filterdir}/cleaned_conserved_overlapping_deletion_conserved_coord_only.txt

#non-reference allele frequency 1, fixed deletions
vcftools --vcf ${vcf_files}/first_test.flt.vcf --keep-only-indels --non-ref-af 1 --recode --out ${vcf_files}/first_test_maf_1 
vcf2bed --deletions < ${vcf_files}/first_test_maf_1.recode.vcf > ${vcf_files}/first_test_maf_1.recode.vcf.bed
#keep only the sites that passed
awk '!($1~"##") && ($8=="PASS") {print}' ${vcf_files}/first_test_maf_1.recode.vcf.bed > ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed

#how many of the fixed deletions overlap with conserved sites?

bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj  | awk '$6!=-1 {print $1"\t"$2"\t"$3}' | sort | uniq | wc -l
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj  | awk '$6!=-1 {print $1"\t"$2"\t"$3}' | sort | uniq | wc -l

#bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj -f 0.9 -r  | awk '$6!=-1 {print}' | head

bedtools closest -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -d  | awk '$NF<10{print}'| wc -l 

#some indels occur away from site? should I locally realign?
#I can take all the reads which annotate the site, and attempt to realign
#what does it mean to have a fixed deletion occuring a bit away from the actual site?
#for conserved sequences overlapping deletions, any deletions from VCF occuring outside of conserved region should be removed
#also need to convert back into original chimp genome locations, as well as human locations

#TODO for conserved overlapping deletions, first see if any complete deletions remove any of the inserted chimp site, then out of these,
#see if these complete deletions overlap the exact annotated site

#converted coordinates first, then original chimp genome coordinates, and human coordinates
awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_hybrid_converted.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_over_only.txt
awk '{print $(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_hybrid_converted.txt | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_over_only.txt

awk '{print $(NF-2)"\t"$(NF-1)"\t"$NF"\t"$0}' ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_hybrid_converted.txt | awk '{for(i=1;i<=NF-3;i++)printf $i"\t"; printf"\n"}' | sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_under_only.txt
awk '{print $(NF-2)"\t"$(NF-1)"\t"$NF"\t"$0}' ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_hybrid_converted.txt | awk '{for(i=1;i<=NF-3;i++)printf $i"\t"; printf"\n"}'| sort -k1,1 -k2,2n | uniq  > ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_under_only.txt


#bedtools closest -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -d  | awk '$NF==0{print}'| wc -l 

#bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_over_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj | awk '$6!=-1 {print}' | cut -f10- | sort -k1,1 -k2,2n | uniq > ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.bed 

######### conservation check ##############

##more broader intersection check##
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_over_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj | awk '$11!=-1 {print}' > ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.txt 

#how many lie completely in conserved region?
awk '$10>$2 && $11<$3 {print $1"\t"$2"\t"$3 }' ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.txt | sort -k1,1 -k2,2n | uniq  | wc -l
#lie partial, its the same as before - good to know
awk '{print $1"\t"$2"\t"$3 }' ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.txt | sort -k1,1 -k2,2n | uniq  | wc -l

##more specific intersection check##
#lie completely in conserved region or not?
#I can probably use this list as a first pass #
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_under_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj | awk '$23!=-1'  > ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.conserved.specific.txt 

#how far are the other fixed deletions from the inserted sites?
#safe bet to simply take fixed deletions lying anywhere in inserted conserved sequence
for i in $(seq 1 30)
do
echo $i
bedtools closest -a ${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_conserved_converted_coord_under_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -d  | awk -v lengthcutoff=$i '$NF<lengthcutoff{print $21"\t"$22"\t"$23}' | sort | uniq | wc -l 
done

##
### for deletions you're only interested if the conserved sequence is completely removed ###, not the other sequences..
###do same thing with deletions###

#do it for the small variants
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_under_only.txt -b ${vcf_files}/first_test_maf_1.flt.vcf.clean.bed -loj | awk '$23!=-1' > ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.deletion.conserved.txt
#deletion overarching count
awk '{print $21"\t"$22"\t"$23}' ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.deletion.conserved.txt | sort | uniq | wc -l
#930

#as well as for the SV 
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_under_only.txt -b ${vcf_files}/first_test.sv.vcf -loj | awk '$23!=-1 && $26=="<DEL>" {print}' > ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.sv.conserved.txt 
awk '{print $21"\t"$22"\t"$23}' ${vcf_files}/first_test_maf_1.flt.vcf.clean.underlie.sv.conserved.txt | sort | uniq | wc -l
#check if there is overlap above?
#576

#for deletions overlapping conserved elements - see if the exact conserved site is deleted

#keep all the VCF files from each genome file id run
#for conserved regions, keep all sites which lie completely in the conserved regions
#for deleted regions, keep all sites in which the deletion completely overlaps the removed conserved regions
#then merge all of them in the end


### structural variants check ###
vcf_files=/cluster_path_temp/sgdp/vcf_files
filterdir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

vcf2bed --deletions < ${vcf_files}/first_test.sv.vcf > ${vcf_files}/first_test.sv.vcf.bed
bedtools intersect -a ${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_conserved_converted_coord_only.txt -b ${vcf_files}/first_test.sv.vcf -loj | awk '$5!=-1{print}' | head

bedtools intersect -a ${chimp_human_hybrid_dir}/chimp_human_hybrid_seq_1_meta.txt -b ${vcf_files}/first_test.sv.vcf -loj | awk '$5!=-1{print}' | head
bedtools intersect -a ${chimp_human_hybrid_dir}/chimp_human_hybrid_seq_1_meta.txt -b ${vcf_files}/first_test.sv.vcf -loj | awk '$8!=-1 && $11=="<DEL>" {print}' |  wc -l

########## bulk analyses ###############

#from the intersected regions, can extract the original human coordinate, original coordinate on chimp genome to do bulk analyses
#bulk analysis thoughts:
#how many of them overlap cortex h3k27ac peaks?#
#how many of them sit close to differentially expressed genes?#
#how many of them overlap 
#for the very small sites, check to see if there is 100% sequence identity?
#this depends on which cell lines are going to be used
#how many of them overlap promoters
#with a list of sites, how do I pick which ones to do for experiments?
#can I test all of them at once? if say I have 10,000 sites?
#how many of them overlap exons?
#do I check to see if they overlap orthologous exons - orthologous as defined by chimp mapping?
#intersect with CDS
awk '{print chr$1"\t"$2"\t"$3}' Chimp_2.1.4_CDS_ensembl_sites.txt | sort -k1,1 -k2,2n | uniq > Chimp_2.1.4_CDS_ensembl_sites.bed
#what cell lines do I use? that determines what functional datasets I eventually overlap

#how many of them overlap defined TF motifs?

#################################################
#play around with vcf files
vcf_files=/cluster_path_temp/sgdp/vcf_files
vcf2bed --deletions < ${vcf_files}/first_test.flt.vcf > ${vcf_files}/first_test.flt.vcf.bed
intersectBed -a -b > 
vcftools  --vcf ${vcf_files}/first_test.flt.vcf --keep-only-indels --out first_test 
vcftools  --vcf ${vcf_files}/first_test.flt.vcf --keep-only-indels --recode --maf 1 --out first_test_maf_1 


vcftools  --vcf ${vcf_files}/first_test.flt.vcf --counts --freq --recode --out ${vcf_files}/first_test_maf_1 


vcftools  --vcf ${vcf_files}/first_test.flt.vcf --keep-only-indels --counts --freq --out first_test_maf_1 
awk '!($1~"##") && ($7=="PASS") {print}' ${vcf_files}/first_test.flt.vcf
#awk -F'[\t:]' 'NR>1{print $1"\t"$2"\t"$6"\t"$8}' first_test_maf_1.frq.count | awk '$4>11{print}'| head

#for indels, if i don't have --keep-only-indels, it only gives me snps?
#awk -F'[\t:]' 'NR>1{print $1"\t"$2"\t"$6"\t"$8}' ${vcf_files}/first_test_maf_1.frq.count | awk '$NF>21{print}' | wc -l
#awk -F'[\t:]' 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8}' ${vcf_files}/first_test_maf_1.frq | awk '($NF==1) && ($4==22){print}' |head

./vcftools --vcf ${vcf_files}/first_test.flt.vcf --non-ref-af 1 --recode --out ${vcf_files}/first_test
awk '!($1~"##") && ($7=="PASS") {print}' ${vcf_files}/first_test.recode.vcf | head
grep chr21 ${vcf_files}/first_test.recode.vcf | grep 46695141

chr21	46695141
awk '$7=="PASS"{print}'  first_test.recode.vcf | wc -l  
#check to see if completely deleted variant can be obtained?
#check to see overlap with conserved sites first? then deleted sites
#####do a few checks to see if converted hybrid genome coordinates match original coordinates

human_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human
filter_dir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

#first three coordinates are new coordinates, other are old
awk '{print $1"\t"$2"\t"$3"\t"$4"|"$5"|"$6}' ${chimp_human_hybrid_dir}/chimp_human_hybrid_seq_1_meta.txt > ${chimp_human_hybrid_dir}/temp_check_ids.txt
bedtools getfasta -name -fi ${chimp_human_hybrid_dir}/chimp_human_hybrid_seq_1.fa -bed ${chimp_human_hybrid_dir}/temp_check_ids.txt -fo ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa

head ${filter_dir}/final_chimp_human_hybrid_seq.fa 
head ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa

diff ${filter_dir}/final_chimp_human_hybrid_seq.fa ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa | head 


#look at chromosomes which differ between two coordinates
awk '$1!=$4{print $1}' ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt | awk -F'[_\t]' '{print $1}' > ${filter_dir}/first_chr_id.txt

awk '$1!=$4{print $4}' ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt | awk -F'[_\t]' '{print $1}' > ${filter_dir}/second_chr_id.txt
paste ${filter_dir}/first_chr_id.txt ${filter_dir}/second_chr_id.txt > ${filter_dir}/all_ids.txt
 awk '$1!=$2{print}' ${filter_dir}/all_ids.txt | uniq -c
 
#look at the first entry to see whats wrong 
#awk 'NR==10655{print}'  ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa
#grep 6468863 ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt
#grep 6257530 ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt

awk 'NR>=10653{print}' ${filter_dir}/final_chimp_human_hybrid_seq.fa | head
awk 'NR>=10653{print}' ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa | head

grep -n "chr11|6237463|6237623" ${filter_dir}/final_chimp_human_hybrid_seq.fa 

#awk '$5=="6468863" {print NR}'  ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt
#awk 'NR==5327{print}'  ${filter_dir}/final_chimp_human_hybrid_replace_coord.txt
grep "chr11|6257530|6257534" ${chimp_human_hybrid_dir}/temp_check_ids.txt

#chr11	6468863	6468867
echo -e "chr11\t6468850\t6468890" > ${chimp_human_hybrid_dir}/temp_check_ids_2.txt
bedtools getfasta -name -fi ${chimp_human_hybrid_dir}/chimp_human_hybrid_seq_1.fa -bed  ${chimp_human_hybrid_dir}/temp_check_ids_2.txt -fo ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa

rm  ${chimp_human_hybrid_dir}/temp_check_ids_2.txt ${chimp_human_hybrid_dir}/temp_check_ids.txt ${chimp_human_hybrid_dir}/final_chimp_human_hybrid_seq_check.fa
