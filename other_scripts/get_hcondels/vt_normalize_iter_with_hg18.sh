gagp_file_header=$1
hcondel_vcf_file_header=$2
gagp_out_path=$3

#for testing
#blacklist_work_dir=/cluster_path/ape_project/deletions_project/blacklist
#gagp_file_header=Pan_troglodytes_FINAL_INDEL
#af_filter=1
#hcondel_vcf_file_header=${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized
#gagp_out_path=/cluster_path_temp/hcondel_blacklist

#consider variant normalization on original file if not variant normalized already?
#may have implications with liftOver, for example if one variant isn't normalized and lies at edge of a sequence

#download first
wget -q -O ${gagp_out_path}/${gagp_file_header}.vcf.gz https://eichlerlab.gs.washington.edu/greatape/data/VCFs/Indels/${gagp_file_header}.vcf.gz 

#normalize, for comparison with hg18-mapped hcondel coordinates later, added 8/29/21 PASS Filter
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${gagp_out_path}/${gagp_file_header}.insertions.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | awk -F"\t" '!($0 ~ /\#/) && (length($5)>length($4)){print}' >> ${gagp_out_path}/${gagp_file_header}.insertions.vcf 

genome_file=/cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.fa
vt normalize ${gagp_out_path}/${gagp_file_header}.insertions.vcf -r ${genome_file} -o ${gagp_out_path}/${gagp_file_header}.hg18.normalized.vcf -w 1000000000  2> ${gagp_out_path}/${gagp_file_header}.hg18.normalized.log
awk -F"\t" '!($0 ~ /\#/){print}' ${gagp_out_path}/${gagp_file_header}.hg18.normalized.vcf > ${gagp_out_path}/${gagp_file_header}.hg18.normalized.header_removed.vcf

#get vcf coordinates information
awk -F"\t" '{print $1"\t"$2-1"\t"$2}' ${gagp_out_path}/${gagp_file_header}.hg18.normalized.header_removed.vcf > ${gagp_out_path}/${gagp_file_header}.vcf.coord_info.txt 

#get allele seq info
awk -F"\t" '{print $4"\t"$5"\t"$7}' ${gagp_out_path}/${gagp_file_header}.hg18.normalized.header_removed.vcf > ${gagp_out_path}/${gagp_file_header}.vcf.allele_seq_info.txt

#get allele frequencies
awk -F"\t" '{print $8}' ${gagp_out_path}/${gagp_file_header}.hg18.normalized.header_removed.vcf | awk -F";" '{print $1"\t"$2"\t"$3}' | sed 's#AC=##g; s#AF=##g; s#AN=##g' > ${gagp_out_path}/${gagp_file_header}.vcf.allele_freq_info.txt 

#rearrange files
paste <(cat ${gagp_out_path}/${gagp_file_header}.vcf.coord_info.txt) <(cat ${gagp_out_path}/${gagp_file_header}.vcf.allele_seq_info.txt) <(cat ${gagp_out_path}/${gagp_file_header}.vcf.allele_freq_info.txt) | awk '{print $1"|"$2"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9"\t0\t+"}' > ${gagp_out_path}/${gagp_file_header}.vcf.allele_info.txt

#added 8/29/21 sort | uniq as vt normalize on GAGP VCF produces some duplicates (with different qualities also it seems)
paste <(cat ${gagp_out_path}/${gagp_file_header}.vcf.coord_info.txt) <(cat ${gagp_out_path}/${gagp_file_header}.vcf.allele_info.txt) | sort | uniq > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt

#added 9/2/21 reset the normalized hg18 VCF for overlap later
awk '{print $4}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | awk -F"|" '{print $1"\t"$2+1"\t"$0"\t"$3"\t"$4"\t0\t"$5"\tNA"}' > ${gagp_out_path}/${gagp_file_header}.hg18.normalized.vcf
 
########## lift to panTro4 (here can choose 2 methods: lift to panTro4, or lift to hg38)  ###############

for genome_ref in panTro4 hg38
do

	if [[ ${genome_ref} == "panTro4" ]]
	then
		liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/hg18ToPanTro4.over.chain.gz
	elif [[ ${genome_ref} == "hg38" ]]
	then
		liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/hg18ToHg38.over.chain.gz
	else
		echo "ref genome must be panTro4 or hg38"
	fi

	genome_file=/cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/${genome_ref}.fa
	chr_size_file=/cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/${genome_ref}.chrom.sizes
    hcondel_vcf_file=${hcondel_vcf_file_header}.${genome_ref}.vcf

	liftOver -minMatch=0.001 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt ${liftChainFile} ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.txt ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.unmapped.txt
	
	#rearrange files
	#https://stackoverflow.com/questions/59472326/what-does-this-mean-awk-ofs-t-1-11-filepath
	#https://stackoverflow.com/questions/16203336/simple-awk-command-issue-fs-ofs-related
	#chr, start pos, strand, followed by meta info
	paste <(awk '{print $1"\t"$2"\t"$6}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.txt) <(awk '{print $4}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.txt) | awk 'BEGIN{FS="|";OFS="\t"} {$1=$1}1' > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.revised.txt 
	
	if [[ ${genome_ref} == "panTro4" ]]
    then 
		#add or subtract ALT sequence length to get insertion annotated
		awk '{if($3=="-"){print $1"\t"$2-length($7)"\t"$2"\t-"} else{ print $1"\t"$2"\t"$2+length($7)"\t+" }}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.revised.txt > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.1.txt
	else
		#for hg38, subtract ALT sequence length if inverted, to get insertion annotated
        awk '{if($3=="-"){print $1"\t"$2-length($7)"\t"$2"\t-"} else{ print $1"\t"$2"\t"$2+1"\t+" }}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.revised.txt > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.1.txt
	fi
	
	paste <(cat ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.1.txt) <(awk '{print $4}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.txt) > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.2.txt 
	#remove sequences less than 0 or greater than chr length
	#NOTE: check discarded sequences later to see if it matches any hcondels to account for edge case of hcondels dropping out due to variant normalization issues
	python /cluster_path/ape_project/deletions_project/check_coord_with_chr_sizes.py ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.2.txt ${chr_size_file} ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.discarded.txt ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.txt

	#remove strand for getfasta
	cut -f1-3,5 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.txt > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.strand_removed.txt

	#get sequences
	bedtools getfasta -name -fi ${genome_file} -bed ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.strand_removed.txt -fo ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa

	#make vcf
	awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa | sed 's#>##g' > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep

	if [[ ${genome_ref} == "panTro4" ]]
    then
		awk -F"[|\t]" '{print $1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"\t"substr($9,1,1)"\t"$9"\t"$5}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep | sort -k1,1 > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep.meta
	
		#make initial VCF, remember that VCF is one based, so add one to coordinate, put strand information at the end
		join -15 -21 <(cat ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.txt | sort -k5,5 ) <(cat ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$6"\t"$7"\tNA\t"$8"\tNA\t"$5}' > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp

		#filter out coordinates whose ${genome_ref}-mapped derived alleles don't match the original read alleles
		#https://stackoverflow.com/questions/62008503/create-reverse-complement-sequence-based-on-awk
		awk -f /cluster_path/ape_project/deletions_project/blacklist/scripts/rev_complement.awk ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.2
		awk -F'[|\t]' 'substr($6,2,length($6))==$17 {print}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.2 > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.filtered

	else
		#hg38, use alt allele from gagp VCF for insertion sequence
		awk -F'[|\t]' '{print $1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"\t"substr($9,1,1)"\t"substr($4,2,length($4))"\t"$5}' ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep | sort -k1,1 > ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep.meta
			
        #make initial VCF, remember that VCF is one based, so add one to coordinate, put strand information at the end
        join -15 -21 <(cat ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.txt | sort -k5,5 ) <(cat ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$6"\t"$7"\tNA\t"$8"\tNA\t"$5}' > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp
			
		#rev comp any alleles that mapped to different strand, and add on bp before	
		awk -f /cluster_path/ape_project/deletions_project/blacklist/scripts/rev_complement_hg38.awk ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.2
		awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$4$10"\t"$6"\t"$7"\t"$8}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.2 > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.filtered
	fi
	
	#make final VCF
	cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf
	cut -f1-8 ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.filtered >> ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf

	#vt normalize to compare
	vt normalize ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf -r ${genome_file} -o ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.vcf -w 1000000000 2> ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalize.log

	#### vt normalize using custom script #####
	#8/9/21 did test and saw it provided same result as vt normalize
	awk '!($0 ~ /\#/) {print $1"\t"$2"\t"$2"\t+\t"$3}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only

	python /cluster_path/ape_project/deletions_project/extendCoordinatesV2_upstream_downstream_vcf.py ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only ${chr_size_file} 3000 100 ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended

	#subtract one for bed format for getfasta
	awk '{print $1"\t"$2-1"\t"$3"\t"$NF}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended.revised
	 
	bedtools getfasta -name -fi ${genome_file} -bed ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended.revised -fo ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq

	awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq | sed 's#>##g' > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised

	#add one back to ref seq coordinate
	paste <( cat ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised ) <( cat ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended.revised) | awk '{print $1"\t"$2"\t"$3"\t"$4+1"\t"$5}' > ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised.with_pos

	#sanity check that paste ids are matched
	#paste <( cat ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised ) <( cat ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended.revised) | awk '$1!=$6{print}'

	python /cluster_path/ape_project/deletions_project/blacklist/scripts/variant_norm_left.py ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised.with_pos ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.custom.vcf 

	#compare to see if everything matches
	vcf_diff=$(diff <(awk '{print $1"_"$2"_"$3"_"$4"_"$5}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.custom.vcf | sort -k1,1) <(awk '!($0 ~ /\#/) {print $1"_"$2"_"$3"_"$4"_"$5}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.vcf | sort -k1,1))
	if [[ ${vcf_diff} != "" ]]
	then
		echo "VCF diff detected"
		echo ${gagp_out_path}/${gagp_file_header}.${genome_ref}
	fi
	
	############################################

	#filter vcf coordinates that pass, join with hcondel vcf
	join -eNA -o auto -t $'\t' -a 1 -11 <( awk '!($0 ~ /\#/) {print $1"_"$2"_"$4"_"$5"\t"$3}' ${hcondel_vcf_file} | sort -k1,1) -21 <( awk '!($0 ~ /\#/) && ($7=="PASS" || $7=="TRF") {print $1"_"$2"_"$4"_"$5"\t"$3}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.vcf | sort -k1,1 ) | awk '{print $2"\t"$3}' > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt

	awk '{print $1}' ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt
	awk '{print $2}' ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt | awk -F"[\t|]" '{if($1=="NA"){print "NA\tNA\tNA"} else{print $6"\t"$7"\t"$8}}' > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt

	paste <(cat ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt) <(cat ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt) > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.txt

done

#mapping from hcondel file to hg18
genome_ref=hg18
hcondel_vcf_file=${hcondel_vcf_file_header}.${genome_ref}.vcf

#filter vcf coordinates that pass, join with hcondel vcf
join -eNA -o auto -t $'\t' -a 1 -11 <( awk '!($0 ~ /\#/) {print $1"_"$2"_"$4"_"$5"\t"$3}' ${hcondel_vcf_file} | sort -k1,1) -21 <( awk '!($0 ~ /\#/) && ($7=="PASS" || $7=="TRF") {print $1"_"$2"_"$4"_"$5"\t"$3}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.vcf | sort -k1,1 ) | awk '{print $2"\t"$3}' > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt

awk '{print $1}' ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt
awk '{print $2}' ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt | awk -F"[\t|]" '{if($1=="NA"){print "NA\tNA\tNA"} else{print $6"\t"$7"\t"$8}}' > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt

paste <(cat ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt) <(cat ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt) > ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.txt

<<COMMENT1

#remove temporary files
rm -f ${gagp_out_path}/${gagp_file_header}.vcf.gz
rm -f ${gagp_out_path}/${gagp_file_header}.hg18.normalized.vcf
rm -f ${gagp_out_path}/${gagp_file_header}.hg18.normalized.header_removed.vcf
rm -f ${gagp_out_path}/${gagp_file_header}.insertions.vcf
rm -f ${gagp_out_path}/${gagp_file_header}.vcf.coord_info.txt
rm -f ${gagp_out_path}/${gagp_file_header}.vcf.allele_seq_info.txt
rm -f ${gagp_out_path}/${gagp_file_header}.vcf.allele_freq_info.txt
rm -f ${gagp_out_path}/${gagp_file_header}.vcf.allele_info.txt

for genome_ref in panTro4 hg38
do
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.mapped.revised.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.1.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.2.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.strand_removed.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep
	rm -f ${gagp_out_path}/${gagp_file_header}.insertion_parsed.${genome_ref}.revised_coord.fa.tab_sep.meta
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.2
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.temp.filtered
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf

	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.coord_only.extended.revised
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised
	rm -f ${gagp_out_path}/${gagp_file_header}.${genome_ref}.vcf.extended_seq.revised.with_pos

	rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt
	rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt
done

genome_ref=hg18
rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.txt	
rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.1.txt
rm -f ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.temp.2.txt

COMMENT1
 
