#8/26/21 checked entire script for typos
#script creating blacklist of all variants

#preprocess
ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.chrom.sizes /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/panTro4.chrom.sizes
ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/panTro4.fa
ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa.fai /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/panTro4.fa.fai

ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.chrom.sizes /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg38.chrom.sizes
ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.fa /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg38.fa
ln -s /cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.fa.fai /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg38.fa.fai

#reset
#blacklist_work_dir=/cluster_path/ape_project/deletions_project/blacklist
#af_filter=1
#rm -f ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered*
#rm -fr /cluster_path_temp/hcondel_blacklist

blacklist_work_dir=/cluster_path/ape_project/deletions_project/blacklist
mkdir -p ${blacklist_work_dir}

gagp_meta_file=${blacklist_work_dir}/gagp/gagp_meta_file.txt

#********************* getting polymorphism data from great ape genome diversity project *********************#

#from GAGP paper
#SNPs were called using GATK after BWA mapping to the human genome (NCBI Build 36) using relaxed mapping parameters.
#Genomes were mapped to the human reference assembly NCBI Build 36 (UCSC hg18) using the BWA mapping software

#download all gagp data first 
mkdir -p ${blacklist_work_dir}/gagp

echo -e "Gorilla_FINAL_INDEL" > ${blacklist_work_dir}/gagp/gagp_meta_file.txt
echo -e "Pan_paniscus_FINAL_INDEL" >> ${blacklist_work_dir}/gagp/gagp_meta_file.txt
echo -e "Pan_troglodytes_FINAL_INDEL" >> ${blacklist_work_dir}/gagp/gagp_meta_file.txt
echo -e "Pongo_abelii_FINAL_INDEL" >> ${blacklist_work_dir}/gagp/gagp_meta_file.txt
echo -e "Pongo_pygmaeus_FINAL_INDEL" >> ${blacklist_work_dir}/gagp/gagp_meta_file.txt

##### normalize hCONDEL file #######
#panTro4 method

#download liftOver
#wget -O /cluster_path/ape_project/deletions_project/liftover_files/hg18ToPanTro4.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToPanTro4.over.chain.gz
#wget -O /cluster_path/ape_project/deletions_project/liftover_files/hg18ToHg38.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz

#get actual bases from fasta

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1
genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa

#all_files_af_${af_filter}_final_filtered.txt is what I used to go into design_oligos.sh
#-1 on chimp del position to get base beforehand for VCF
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt

paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2-1"\t"$3"\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.txt

#get it in VCF format, then vt normalize

#get sequences
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa

#make vcf
awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa.tab_sep

awk '{print $1"\t"substr($2,1,1)"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa.tab_sep | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa.tab_sep.meta

#vcf one based, so add one to start coordinate
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.panTro4.vcf
join -14 -21 <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.txt | sort -k4,4 ) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.panTro4.chimp_del.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$5"\t"$6"\tNA\tPASS\tNA"}' >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.panTro4.vcf

vt normalize ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.panTro4.vcf -r ${genome_file} -o ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf -w 1000000000  2> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.log


##### normalize hCONDEL file #######
#hg38 method

#download liftOver
#wget -O /cluster_path/ape_project/deletions_project/liftover_files/hg18ToPanTro4.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToPanTro4.over.chain.gz

#get actual bases from fasta

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1

#all_files_af_${af_filter}_final_filtered.txt is what I used to go into design_oligos.sh
#get the panTro4 sequence for the deletion, and use the hg38 base before the insertion
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt
paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.txt

#hg38 original coordinate in 1-based, so subtract one to get bed
paste <(cut -f10-11 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.txt

#get it in VCF format, then vt normalize

#get sequences
genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa

genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.fa
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.fa

#make vcf
#get strand
awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa | sed 's#>##g' | awk -F'[\t|]' '{print $0"\t"$13}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep

#reverse complement
awk -f /cluster_path/ape_project/deletions_project/blacklist/scripts/rev_complement_v2.awk ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep >  ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp

awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.fa | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.fa.tab_sep

#concatenate human base before insertion and chimp insertion sequences 
paste <(awk '{print $1"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.fa.tab_sep) <(awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp) | awk '{print $1"\t"$2"\t"$2$3}' | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_human_del.fa.tab_sep.meta

#paste <(awk '{print $1}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.fa.tab_sep) <(awk '{print $1}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp) | awk '$1!=$2 {print}'
#awk '{print $1"\t"substr($2,1,1)"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.meta

#vcf one based, so add one to coordinate
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg38.vcf
join -14 -21 <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.txt | sort -k4,4 ) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_human_del.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$5"\t"$6"\tNA\tPASS\tNA"}' >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg38.vcf

genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.fa
vt normalize ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg38.vcf -r ${genome_file} -o ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf -w 1000000000  2> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.log

########## map hcondel file to hg18 also as a sanity check ##########
#already ran before
#paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.txt

#hg38 original coordinate in 1-based, so subtract one to get bed
paste <(cut -f10-11 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.txt

#map to hg19
#wget -O /cluster_path/ape_project/deletions_project/liftover_files/hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg38ToHg19.over.chain.gz
liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/hg38ToHg19.over.chain.gz
liftOver -minMatch=0.001 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg19_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg19_unmapped.txt

#map to hg18
#wget -O /cluster_path/ape_project/deletions_project/liftover_files/hg19ToHg18.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg19ToHg18.over.chain.gz
liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/hg19ToHg18.over.chain.gz
liftOver -minMatch=0.001 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg19_mapped.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_unmapped.txt

#already ran before, get sequences
#genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa
#bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa

genome_file=/cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.fa
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa

#make vcf
#already ran before, get strand
#awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa | sed 's#>##g' | awk -F'[\t|]' '{print $0"\t"$13}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep

#already ran before, reverse complement
#awk -f /cluster_path/ape_project/deletions_project/blacklist/scripts/rev_complement_v2.awk ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep >  ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp

awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa.tab_sep

#concatenate human base before insertion and chimp insertion sequences 
#old
#paste <(awk '{print $1"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa.tab_sep) <(awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp) | awk '{print $1"\t"$2"\t"$2$3}' | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.chimp_human_del.fa.tab_sep.meta
join -11 -21 <(awk '{print $1"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa.tab_sep | sort -k1,1 ) <(awk '{print $1"\t"$4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp | sort -k1,1) | awk '{print $1"\t"$2"\t"$2$3}' | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.chimp_human_del.fa.tab_sep.meta

#paste <(awk '{print $1}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.human_del.fa.tab_sep) <(awk '{print $1}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.rev_comp) | awk '$1!=$2 {print}'
#awk '{print $1"\t"substr($2,1,1)"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep | sort -k1,1 > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.chimp_del.fa.tab_sep.meta

#vcf one based, so add one to coordinate
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg18.vcf
join -14 -21 <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt | sort -k4,4 ) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg18.chimp_human_del.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$5"\t"$6"\tNA\tPASS\tNA"}' >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg18.vcf

#wget -q -O /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit 
#twoBitToFa /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.2bit /cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.fa

genome_file=/cluster_path/ape_project/deletions_project/blacklist/reference_genome_files/hg18.fa
vt normalize ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.initial.hg18.vcf -r ${genome_file} -o ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg18.vcf -w 1000000000  2> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg18.log


###### submit to cluster to run ##########
hcondel_vcf_file_header=${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized
gagp_out_path=/cluster_path_temp/hcondel_blacklist
mkdir -p ${gagp_out_path}
header_name=gagp_hcondel_process

bash ${blacklist_work_dir}/scripts/vt_normalize_submit_with_hg18.sh ${gagp_meta_file} ${hcondel_vcf_file_header} ${gagp_out_path} ${header_name}

#tail ${gagp_out_path}/${header_name}_vt_normalize_submit_with_hg18.log
#********************* mapping to other reference genomes  *********************#
#DOUBLE CHECK
#map to panTro5, see how many are deleted
#left align, right align then map

#compare sequences where vcf annotates it as a complete deletion but not pantro5

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1
genome_file=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.fa
chr_size_file=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp/panTro4.chrom.sizes
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt
#all_files_af_${af_filter}_final_filtered.txt is what I used to go into design_oligos.sh

#get extended seq positions for left/right alignment
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt
paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2"\t"$3"\t+\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt

python /cluster_path/ape_project/deletions_project/extendCoordinatesV2_upstream_downstream.py ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt ${chr_size_file} 1000 1000 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt.extended

awk '{print $1"\t"$2"\t"$3"\t"$NF}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt.extended > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt.extended.revised

bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt.extended.revised -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended

awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended.revised

#add one back to ref seq coordinate
paste <( cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended.revised ) <( cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.txt.extended.revised) | awk '{print $1"\t"$2"\t"$3"\t"$4+1"\t"$5}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended.revised.with_pos


################################# left align #################################
#subtract one to get sequence before
paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2-1"\t"$3"\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.txt

#get it in VCF format, then vt normalize

#get sequences
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa

#make vcf
awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa.tab_sep

awk '{print $1"\t"substr($2,1,1)"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa.tab_sep > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa.tab_sep.meta

#vcf one based, so add one to coordinate
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.initial.vcf 
paste <(awk '{print $4"\t"$1"\t"$2"\t"$3}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.txt ) <(cut -f2- ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.fa.tab_sep.meta) | awk '{print $2"\t"$3+1"\t"$1"\t"$5"\t"$6"\tNA\tPASS\tNA"}' >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.initial.vcf

python /cluster_path/ape_project/deletions_project/blacklist/scripts/variant_norm_left.py ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.initial.vcf ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended.revised.with_pos ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.vcf

################################# right align #################################

paste <(cut -f7-9 ${deletions_data_coord} ) <( awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} ) | awk '{print $1"\t"$2"\t"$3+1"\t"$4}' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.txt

#get it in VCF format, then vt normalize

#get sequences
bedtools getfasta -name -fi ${genome_file} -bed ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.txt -fo ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa

#make vcf
awk '!/^>/ { next } { getline seq } { print $0"\t"toupper(seq) }' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa | sed 's#>##g' > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa.tab_sep

awk '{print $1"\t"substr($2,length($2),length($2))"\t"$2}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa.tab_sep > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa.tab_sep.meta

#vcf one based, $4-1 is one-base position of variant sequence
cat /cluster_path/ape_project/deletions_project/blacklist/scripts/vcf_file_header.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.initial.vcf 
paste <(awk '{print $4"\t"$1"\t"$2"\t"$3}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.txt ) <(cut -f2- ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.fa.tab_sep.meta) | awk '{print $2"\t"$4"\t"$1"\t"$5"\t"$6"\tNA\tPASS\tNA"}' >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.initial.vcf

python /cluster_path/ape_project/deletions_project/blacklist/scripts/variant_norm_right.py ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.initial.vcf ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.extended.revised.with_pos ${chr_size_file} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.vcf

################################# lift to other reference genomes ################################# 

#wget -O /cluster_path/ape_project/deletions_project/liftover_files/panTro4ToPanTro5.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/panTro4/liftOver/panTro4ToPanTro5.over.chain.gz

#left coord
#$2 is one bp before deleted sequence (one-based)
awk '{print $1"\t"($2+1)-1"\t"($2+1)-1+(length($5)-1)"\t"$3}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.vcf > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.txt

#right coord
#$2 is one bp after deleted sequence (one-based)
awk '{print $1"\t"($2-1)-(length($5)-1)"\t"$2-1"\t"$3}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.vcf > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.txt

while read -a array
do
	genome_ref="${array[0]}"
	genome_ref_cap=$(echo ${genome_ref} | awk '{ print toupper(substr($0, 1, 1)) substr($0, 2) }')
	liftChainFile="${array[1]}"

	#left, lift to ${genome_ref}

	liftOver -minMatch=1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.txt

	grep -v "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.cleaned.txt

	#right, lift to ${genome_ref}

	liftOver -minMatch=1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.txt

	grep -v "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.cleaned.txt

	##########

	#overlap with hcondels deleted from panTro4
	awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.cleaned.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt
	awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.cleaned.txt >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt
	sort -k1,1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt | uniq > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.uniq.txt

	join -e1 -o auto -t $'\t' -a 1 -11 -21 <(awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} | sort -k1,1) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.uniq.txt | awk '{print $1"\t"0}') > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_${genome_ref}_map_blacklist_info.txt

done</cluster_path/ape_project/deletions_project/blacklist/lift_info.txt

################ panTro6 check for reviewers, added 4/20/22 #######################

post_vcf_filter_dir=/cluster_path/ape_project/deletions_project/post_vcf_filtering
af_filter=1
deletions_data_coord=${post_vcf_filter_dir}/all_files_af_${af_filter}_final_filtered.txt

genome_ref=panTro6
genome_ref_cap=$(echo ${genome_ref} | awk '{ print toupper(substr($0, 1, 1)) substr($0, 2) }')
liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/panTro5ToPanTro6.over.chain.gz

#left, lift to ${genome_ref}

liftOver -minMatch=1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.panTro5_mapped.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.txt

grep -v "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.cleaned.txt

#right, lift to ${genome_ref}

liftOver -minMatch=1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.panTro5_mapped.txt ${liftChainFile} ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_mapped.txt ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.txt

grep -v "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.cleaned.txt

#look at unmapped sequences

#grep "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.txt | sort | uniq -c
#grep "#" ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.txt | sort | uniq -c
##########

#overlap with hcondels deleted from panTro4
awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.${genome_ref}_unmapped.cleaned.txt > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt
awk '{print $4}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.${genome_ref}_unmapped.cleaned.txt >> ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt
sort -k1,1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.txt | uniq > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.uniq.txt

join -e1 -o auto -t $'\t' -a 1 -11 -21 <(awk 'BEGIN{OFS="|";} {$1=$1}1' ${deletions_data_coord} | sort -k1,1) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.${genome_ref}_unmapped.cleaned.uniq.txt | awk '{print $1"\t"0}') > ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_${genome_ref}_map_blacklist_info.txt

#awk '$2==0 {print}' ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_${genome_ref}_map_blacklist_info.txt | wc -l

#********************* merge files  *********************#

headerString="SeqName"

for genome_ref in panTro4 hg38 hg18
do
	gagp_file_header=$(awk 'NR==1{print $1}' ${blacklist_work_dir}/gagp/gagp_meta_file.txt)
	gagp_merged_file=${gagp_out_path}/gagp_merged.hcondel_overlapped.${genome_ref}.txt
	cat ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.txt > ${gagp_merged_file}
	
	awk 'NR>1{print}' ${blacklist_work_dir}/gagp/gagp_meta_file.txt > ${blacklist_work_dir}/gagp/gagp_meta_file.revised.txt
	
	headerStringIter="${gagp_file_header}_AC_${genome_ref}\t${gagp_file_header}_AF_${genome_ref}\t${gagp_file_header}_AN_${genome_ref}"
	headerString="${headerString}\t${headerStringIter}"
	
	while read -a array
	do
		gagp_file_header="${array[0]}"
		gagp_iter_file=${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.txt
		wc -l ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.${genome_ref}.txt
		join -t $'\t' -11 -21 <(cat ${gagp_merged_file} | sort -k1,1) <(cat ${gagp_iter_file} | sort -k1,1) >> ${gagp_merged_file}.temp
		mv ${gagp_merged_file}.temp ${gagp_merged_file}
	
		headerStringIter="${gagp_file_header}_AC_${genome_ref}\t${gagp_file_header}_AF_${genome_ref}\t${gagp_file_header}_AN_${genome_ref}"
		headerString="${headerString}\t${headerStringIter}"
		
	done<${blacklist_work_dir}/gagp/gagp_meta_file.revised.txt
done 

gagp_merged_file_final=${blacklist_work_dir}/novaseq_12_15_18_blacklist_info.txt

#merge panTro4, hg38, and hg18 information
join -t $'\t' -11 -21 <(cat ${gagp_out_path}/gagp_merged.hcondel_overlapped.panTro4.txt | sort -k1,1) <(cat ${gagp_out_path}/gagp_merged.hcondel_overlapped.hg38.txt | sort -k1,1) > ${gagp_merged_file_final}.temp
mv ${gagp_merged_file_final}.temp ${gagp_merged_file_final}
join -eNA -o auto -t $'\t' -a 1 -11 <(cat ${gagp_merged_file_final} | sort -k1,1) -21 <(cat ${gagp_out_path}/gagp_merged.hcondel_overlapped.hg18.txt | sort -k1,1) > ${gagp_merged_file_final}.temp
mv ${gagp_merged_file_final}.temp ${gagp_merged_file_final}

#add on other ref genome map annotation
while read -a array
do
	genome_ref="${array[0]}"
	join -t $'\t' -11 -21 <(sort -k1,1 ${gagp_merged_file_final}) <(sort -k1,1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_${genome_ref}_map_blacklist_info.txt) > ${gagp_merged_file_final}.temp
	mv ${gagp_merged_file_final}.temp ${gagp_merged_file_final}
	headerString="${headerString}\tpanTro4_${genome_ref}"
done</cluster_path/ape_project/deletions_project/blacklist/lift_info.txt

#add on panTro6
genome_ref=panTro6
join -t $'\t' -11 -21 <(sort -k1,1 ${gagp_merged_file_final}) <(sort -k1,1 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_${genome_ref}_map_blacklist_info.txt) > ${gagp_merged_file_final}.temp
mv ${gagp_merged_file_final}.temp ${gagp_merged_file_final}
headerString="${headerString}\tpanTro4_${genome_ref}"

#add on header
echo -e "${headerString}" > ${gagp_merged_file_final}.temp
cat "${gagp_merged_file_final}" >> ${gagp_merged_file_final}.temp

mv ${gagp_merged_file_final}.temp ${gagp_merged_file_final}

#contains the above commands
#bash ${blacklist_work_dir}/scripts/merge_gagp_meta_files_hg18_added_panTro6.sh

#********************* now get initial statistics *********************#

#for VCF, lift hg18 to hg38?, variant normalize again?  then compare?

#how many hcondels are deleted in chimp?
awk  '$9==1{print}' ${gagp_merged_file}  | wc -l

#how many hcondels are deleted in chimp and bonobo?
awk  '$6==1 && $9==1{print}' ${gagp_merged_file}  | wc -l

#how many hcondels are deleted in chimp and bonobo and don't have panTro5 issue?
awk  '$6==1 && $9==1 && $NF==1{print}' ${gagp_merged_file}  | wc -l

#how many hcondels are polymorphic in chimp?
awk  '$9!="NA" && $9!=1{print}' ${gagp_merged_file}  | wc -l

#how many hcondels are deleted in all primates?
awk  '$6==1 && $9==1 && $12==1 && $15==1{print}' ${gagp_merged_file}  | wc -l

#17749
#13402 completely deleted
#2445 polymorphic
#1149 fake hcondels

#see how many hcondels are blacklisted

#why are multiple variants in the same place? double check original variants
#check to see if the aggregated chimps add up to all 50 chimps? could be due to varying read depth that certain hcondels captured by all chimps/certain not?

gagp_out_path=/cluster_path_temp/hcondel_blacklist

#which ones are in panTro4, but not in hg38
join -v1 -11 <(awk '$9==1{print $1}' ${gagp_out_path}/gagp_merged.hcondel_overlapped.panTro4.txt | sort -k1,1) <(awk '$9==1{print $1}' ${gagp_out_path}/gagp_merged.hcondel_overlapped.hg38.txt | sort -k1,1) | wc -l

#which ones are in hg38, but not in panTro4
join -v2 -11 <(awk '$9==1{print $1}' ${gagp_out_path}/gagp_merged.hcondel_overlapped.panTro4.txt | sort -k1,1) <(awk '$9==1{print $1}' ${gagp_out_path}/gagp_merged.hcondel_overlapped.hg38.txt | sort -k1,1) | wc -l

#add a column for rhemac8?

#********************* additional checks *********************#

#not found in chimp, not found in pantro4 to pantro5 - fake hCONDELs
join -11 -21 <(awk '$9=="NA" && $9!=1{print $1}'  ${gagp_merged_file} | sort -k1,1) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.panTro5_unmapped.cleaned.uniq.txt | sort -k1,1) | wc -l

#check which sequences are deleted in both left align and right alignment? - these could be gaps in alignment
join -14 -24 <(sort -k4,4 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.left_align.del_coord.panTro5_unmapped.cleaned.txt) <(sort -k4,4 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.right_align.del_coord.panTro5_unmapped.cleaned.txt) | wc -l

#NA not captured by panTro4 to panTro5 map
join -v 1 -11 -21 <(awk '$9=="NA"{print $1}'  ${gagp_merged_file} | sort -k1,1) <(cat ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered_chimp_del.panTro5_unmapped.cleaned.uniq.txt | sort -k1,1) | wc -l

gagp_out_path=/cluster_path_temp/hcondel_blacklist
gagp_file_header=Pan_troglodytes_FINAL_INDEL
af_filter=1

#sequences to lookup in MSA
#chr10|124662962|124662979|chr10|125100702|125100714|chr10|124662964|124662969|chr10|125100704|125100704|+|C
#chr3|9340915|9340943|chr3|9143286|9143302|chr3|9340919|9340931|chr3|9143290|9143290|+|C
#chr10|124662962|124662979|chr10|125100702|125100714|chr10|124662964|124662969|chr10|125100704|125100704|+|C

#left aligned in 30 mammals, but right aligned in 100 vertebrates
#chr5|104364774|104364851|chr5|104004232|104004308|chr5|104364803|104364804|chr5|104004261|104004261|+|C

#chimp-specific insertion - may be false since gagp didn't detect it?
#chr10|100474077|100474336|chr10|101284095|101284353|chr10|100474167|100474168|chr10|101284185|101284185|+|C

#hard to sequence region?, lots of gaps
#chr10|112206445|112207267|chr10|112818406|112819228|chr10|112206501|112206502|chr10|112818462|112818462|+|C
#chr3|83321857|83322003|chr3|81808838|81808982|chr3|83321997|83321999|chr3|81808978|81808978|+|C
#chr3|54532073|54532263|chr3|53665637|53665825|chr3|54532084|54532086|chr3|53665648|53665648|+|C
#chr7|131223424|131223497|chr7|129744398|129744453|chr7|131223476|131223494|chr7|129744450|129744450|+|C

#lookup vcf 
#chr5|181203199|181203235|chr5|180051843|180051872|chr5|181203206|181203213|chr5|180051850|180051850|+|C

#lookup vcf, no idea why not found
#chr17|78757895|78758023|chr17|79593373|79593501|chr17|78758022|78758023|chr17|79593501|79593501|+|C

#panTro4 has a CT being deleted, panTro5 has a CT, VCF has AT as deletion in all primates except panTro4 - thats why it got filtered out
#chr15|35792550|35792671|chr15|38919548|38919667|chr15|35792650|35792652|chr15|38919648|38919648|+|C
grep 357926 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.normalized.vcf | grep chr15

#panTro4 has a T being deleted, panTro5 has a G, VCF has G as deletion, so most likely error in panTro4 alignment?
#chr3|41842087|41842841|chr3|41194787|41195540|chr3|41842301|41842302|chr3|41195001|41195001|+|C
grep 41842 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.normalized.vcf | grep chr3

#panTro4 has a ATAG being deleted, panTro5 has ATAG, but in VCF, possible that deletion is split into multiple pieces (so polymorphic)
#double check later if interested
#chr13|53468552|53468601|chr13|53800674|53800719|chr13|53468594|53468598|chr13|53800716|53800716|+|C
grep 53468 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.normalized.vcf | grep chr13

#panTro4 has two Aa being deleted, panTro 5 has one A being deleted
#chr16|73632228|73632445|chr16|74265321|74265536|chr16|73632233|73632235|chr16|74265326|74265326|+|C

#panTro4 has two C's being deleted, panTro 5 has 2 Cs being deleted, VCF has 1 C being deleted
#chr20|47106498|47106778|chr20|50054973|50055251|chr20|47106666|47106668|chr20|50055141|50055141|+|C
grep 471066 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.normalized.vcf | grep chr20

#panTro4 has 5 G's being deleted, panTro5 has 3 G's
#chr22|25827553|25827710|chr22|27158101|27158253|chr22|25827698|25827703|chr22|27158246|27158246|+|C
#panTro4 has 1 A being deleted, panTro5 has two A's being deleted, VCF not annotated
grep 25827 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.normalized.vcf | grep chr22
grep 25827 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2 | grep chr22

#panTro4 has a A panTro 5 doesn't?, liftover maps A to G
#chr4	42135762	42135763 is deleted, so first G is deleted instead
#chr4|42135672|42135788|chr4|42067082|42067197|chr4|42135765|42135766|chr4|42067175|42067175|+|C

#panTro4 has two Ts, panTro5 has 3 Ts
#chr13|53981195|53981226|chr13|54303529|54303558|chr13|53981221|53981223|chr13|54303555|54303555|+|C

#panTro4 right aligned, panTro5 left aligned on human alignment, can't find in VCF
#chr11|130261450|130261761|chr11|132199655|132199965|chr11|130261542|130261543|chr11|132199747|132199747|+|C
grep 130261 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2  | grep chr11

#panTro4 ccaag deleted, panTro5 gacaa deleted, difference in C to A mutation?, VCF appears to agree with panTro5
#panTro4  tttgtaagccaagacaaaat
#panTro5 tttgtaagacaagacaaaat
#chr9|74398136|74398212|chr9|75799303|75799374|chr9|74398170|74398175|chr9|75799337|75799337|+|C
grep 743981 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2  | grep chr9

#panTro4 right aligned, panTro5 left aligned on human alignment, VCF 2 Ts deleted
#chr14|28241558|28241644|chr14|29382030|29382115|chr14|28241584|28241585|chr14|29382056|29382056|+|C
grep 282415 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2  | grep chr14

#deletion of one A in panTro4 and panTro5, with T snp next to it, can't find in VCF
#chr18|6216780|6216936|chr18|10366152|10366306|chr18|6216789|6216790|chr18|10366297|10366297|-|C
grep 62167 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2  | grep chr18

#panTro4 ccttccagc is deleted, panTro5 gcccttcca is deleted, VCF has GCCCTTCCA, but there is a N in the read which makes it not align (N is C)
#panTro4 ccctgcccttccagctgg------------g
#panTro5 ccctgcccttccagctgg------------g
#chr19|46003786|46003851|chr19|40796535|40796597|chr19|46003848|46003857|chr19|40796597|40796597|+|P
grep 460038 ${gagp_out_path}/${gagp_file_header}.panTro4_mapped.vcf.temp.2  | grep chr19


########### which ones are in panTro4, but not in hg38 ###########

#two different alignments, can have Ts deleted at the end or in the beginning
#chr10|111933190|111933535|chr10|112550229|112550569|chr10|111933264|111933267|chr10|112550303|112550303|+|C

#find the hg18 aligned reads
grep 111933264 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.panTro4.temp.txt
grep 111933264 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 114300045 ${gagp_out_path}/${gagp_file_header}.hg38.vcf
grep 114300045 ${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf


#two different alignments, can go either way
#chr3|51294145|51294264|chr3|50497834|50497951|chr3|51294253|51294255|chr3|50497942|50497942|+|C

#find the hg18 aligned reads
grep 51294253 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.panTro4.temp.txt
grep 51294253 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 51294253 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 50510362 ${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf

#chr8|33112512|33112633|chr8|36631464|36631583|chr8|33112597|33112598|chr8|36631549|36631549|+|C


#two different alignments, can go either way
#chr8|33112512|33112633|chr8|36631464|36631583|chr8|33112597|33112598|chr8|36631549|36631549|+|C

#find the hg18 aligned reads
grep 33112597 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.panTro4.temp.txt
grep 33112597 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 33112597 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 36608218 ${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf


#two different alignments, can go either way
#chr5|21499144|21499765|chr5|93875912|93876529|chr5|21499156|21499160|chr5|93876517|93876517|-|C

#find the hg18 aligned reads
grep 21499156 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.panTro4.temp.txt
grep 21499156 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 21499156 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 93237974 ${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf


#two different alignments, can go either way
#chr4|50512534|50512891|chr4|79526755|79527112|chr4|50512882|50512883|chr4|79526763|79526763|-|C

#find the hg18 aligned reads
grep 50512882 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.panTro4.temp.txt
grep 50512882 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 50512882 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 80666930 ${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf


########### which ones are in hg38, but not in panTro4 ###########

#unmapped from liftOver, very gappy region
#chr10|72233258|72233284|chr10|73635844|73635869|chr10|72233260|72233261|chr10|73635846|73635846|+|C

#find the hg18 aligned reads
grep 72233260 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.hg38.temp.txt
grep 72233260 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 72233260 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 75065608 ${gagp_out_path}/${gagp_file_header}.panTro4.normalized.vcf

grep 75065608 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.panTro4.mapped.revised.txt

grep 75065608 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.panTro4.unmapped.txt

#maps to completely different region from hg18 to panTro4
#chr1|117862231|117862335|chr1|118849246|118849348|chr1|117862296|117862298|chr1|118849283|118849283|-|C

#find the hg18 aligned reads
grep 117862296 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.hg38.temp.txt
grep 117862296 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 117862296 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 119193428 ${gagp_out_path}/${gagp_file_header}.panTro4.normalized.vcf

#maps to completely different region (chr20) from hg18 to panTro4
#chr7|61549499|61549604|chr7|57594220|57594322|chr7|61549545|61549548|chr7|57594266|57594266|+|C

#find the hg18 aligned reads
grep 61549545 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.hg38.temp.txt
grep 61549545 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 61549545 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 57657912 ${gagp_out_path}/${gagp_file_header}.panTro4.normalized.vcf

grep 57657912 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.panTro4.mapped.txt

#unmapped from liftOver, very gappy region nearby
#chr8|32050961|32051118|chr8|35569358|35569489|chr8|32051026|32051028|chr8|35569398|35569398|+|C

#find the hg18 aligned reads
grep 32051026 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.hg38.temp.txt
grep 32051026 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 32051026 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 35546457 ${gagp_out_path}/${gagp_file_header}.panTro4.normalized.vcf

grep 35546457 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.panTro4.mapped.revised.txt

grep 35546457 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.panTro4.unmapped.txt

#maps to completely different region from hg18 to panTro4
#chr6_GL390389_random|601913|602071|chr6|60390771|60390919|chr6_GL390389_random|602046|602053|chr6|60390901|60390901|+|C

#find the hg18 aligned reads
grep 602046 ${gagp_out_path}/${gagp_file_header}.hcondel_overlapped.hg38.temp.txt
grep 602046 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf

grep 602046 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.panTro4.vcf

grep 57470896 ${gagp_out_path}/${gagp_file_header}.panTro4.normalized.vcf

########### which ones are not in hg38 or panTro4, but in panTro5 ###########
#chr1|104426711|104427022|chr1|103415381|103415689|chr1|104426996|104426999|chr1|103415666|103415666|+|C
grep 10341 ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt

########### which ones are not in VCF and not missing in panTro5 ###########
#these are essentially indels without VCF support

#chr1|104426711|104427022|chr1|103415381|103415689|chr1|104426996|104426999|chr1|103415666|103415666|+|C
#not present
grep 103415666 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 103415666 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b10373" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr1\b"

#chr3|115240490|115240541|chr3|111924471|111924517|chr3|115240510|115240515|chr3|111924491|111924491|+|C
#GCTTC (VCF) vs GCTAC (original) disagreement
grep 111924491 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 111924491 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b11312" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr3\b"

#chr10|100474077|100474336|chr10|101284095|101284353|chr10|100474167|100474168|chr10|101284185|101284185|+|C
#not present
grep 101284185 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 101284185 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b10303" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr10\b"

#chr11|100450630|100450658|chr11|102606022|102606042|chr11|100450650|100450658|chr11|102606042|102606042|+|C
#difference in alleles, CATC (VCF) vs CATCCATC (original)
grep 102606042 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 102606042 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b10198" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr11\b"

#chr19|10649158|10649203|chr19|10454518|10454544|chr19|10649179|10649198|chr19|10454539|10454539|+|C
#difference in alleles, TTCCTCAGGT (VCF) vs GCCCCTCCCTTCCACAGGT (original)
grep 10454539 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 10454539 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b10426" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr19\b"

#chr20|7327894|7327979|chr20|7437453|7437534|chr20|7327972|7327976|chr20|7437531|7437531|+|C
#not present, but nearby has split polymorphism 7366180, A, 7366182, TGT
grep 7437531 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 7437531 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b7366" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr20\b"

#chr22|33477206|33477233|chr22|34795585|34795610|chr22|33477209|33477211|chr22|34795588|34795588|+|C
#difference in alleles, GA (VCF) vs GAA (original)
grep 34795588 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 34795588 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b33521" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr22\b"

#chr7|134660793|134660824|chr7|133165653|133165682|chr7|134660799|134660801|chr7|133165659|133165659|+|C
#difference in alleles, CTT (VCF) vs CT (original)
grep 133165659 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 133165659 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b13250" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr7\b"

#chr15|91179985|91180142|chr15|93775012|93775165|chr15|91180133|91180137|chr15|93775160|93775160|+|C
#difference in alleles, CA (VCF) vs CAAAA (original), 90% frequency also
grep 93775160 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 93775160 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b92119" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr15\b"

#chr17|38269602|38269800|chr17|17681270|17681463|chr17|38269605|38269606|chr17|17681460|17681460|-|C
#not present
grep 17681460 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 17681460 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b17525" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr17\b"

#positive control 1
#chr12|33386604|33386700|chr12|55977690|55977779|chr12|33386649|33386656|chr12|55977734|55977734|-|C 
grep 55977734 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 55977734 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b54657" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr12\b"

#positive control 2
#chr4|26290856|26290974|chr4|26321776|26321891|chr4|26290953|26290956|chr4|26321873|26321873|+|C 
grep 26321873 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.hg38.human_del.hg18_mapped.txt  
grep 26321873 ${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized.hg38.vcf
grep "\b25932" ${gagp_out_path}/${gagp_file_header}.insertion_parsed.txt | grep "\bchr4\b"

${gagp_out_path}/${gagp_file_header}.hg38.normalized.vcf

########### which ones are in new blacklist, but not in old blacklist? ###########
#not normalized, no PASS filter used, some duplicates in there
#chr10|35373783|35373904|chr10|34896723|34896840|chr10|35373870|35373874|chr10|34896810|34896810|+|C
grep 35373874 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 35225744

#chr13|105175767|105175856|chr13|104858960|104859045|chr13|105175793|105175797|chr13|104858986|104858986|+|C
grep 105175797 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 104309338

#chr16|62317251|62317512|chr16|63133195|63133454|chr16|62317372|62317373|chr16|63133316|63133316|+|C
grep 62317373 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 61724715

#chr17|8982389|8982578|chr17|48615315|48615503|chr17|8982516|8982517|chr17|48615376|48615376|-|C
grep 8982517 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 44047732

#chr1|149154520|149154653|chr1|170737085|170737216|chr1|149154617|149154619|chr1|170737182|170737182|+|C
grep 149154619 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 168972934

#chr9|122872481|122872608|chr9|123814022|123814146|chr9|122872534|122872537|chr9|123814075|123814075|+|C
grep 122872537 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 125616170

#chr6|35579390|35579638|chr6|35137352|35137598|chr6|35579508|35579510|chr6|35137470|35137470|+|C
grep 35579510 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 35213218

#chr4|14011902|14012083|chr4|14190260|14190436|chr4|14011954|14011955|chr4|14190312|14190312|+|C
grep 14011955 ${hcondel_vcf_file_header}.${genome_ref}.vcf
gunzip -c ${gagp_out_path}/${gagp_file_header}.vcf.gz | grep 14011955


########### which ones are in new blacklist, but not in old blacklist? ###########
blacklist_work_dir=/cluster_path/ape_project/deletions_project/blacklist
genome_ref=panTro4
gagp_out_path=/cluster_path_temp/hcondel_blacklist
gagp_file_header=Pan_troglodytes_FINAL_INDEL
af_filter=1
hcondel_vcf_file_header=${blacklist_work_dir}/gagp/all_files_af_${af_filter}_final_filtered.normalized

#chr3|17962800|17963122|chr3|17706872|17707187|chr3|17962804|17962806|chr3|17706876|17706876|+|C 
#QDFilter in old removed example

grep 17962806 ${hcondel_vcf_file_header}.${genome_ref}.vcf
grep 1772336 ${gagp_out_path}/${gagp_file_header}.insertions.vcf | grep chr3

#chr3|186953047|186953132|chr3|182943500|182943581|chr3|186953047|186953051|chr3|182943500|182943500|+|C 

genome_ref=panTro4
grep 186953047 ${hcondel_vcf_file_header}.${genome_ref}.vcf

genome_ref=hg18
grep 186953047 ${hcondel_vcf_file_header}.${genome_ref}.vcf

grep 184143978 ${gagp_out_path}/${gagp_file_header}.insertions.vcf 

#chr10|19528858|19529002|chr10|19280109|19280252|chr10|19528915|19528916|chr10|19280166|19280166|+|C 

genome_ref=panTro4
grep 19528916 ${hcondel_vcf_file_header}.${genome_ref}.vcf

genome_ref=hg18
grep 19528916 ${hcondel_vcf_file_header}.${genome_ref}.vcf

grep 19609095 ${gagp_out_path}/${gagp_file_header}.insertions.vcf 


#####

#check that all ids are different are due to a filter being off
while read -a array
do
	gagp_file_header="${array[0]}"
	echo "${gagp_file_header}"
	
	#get hg18 coordinates
	genome_ref=hg18
	join -eNA -o auto -t $'\t' -a 1 -11 <(awk '{print}' ${blacklist_work_dir}/discord_files/novaseq_12_15_18_${gagp_file_header}_discordant_ids.txt | sort -k1,1) -21 <(awk '!($0 ~ /\#/) {print $3"\t"$1"_"$2"_"$4"_"$5}' ${hcondel_vcf_file_header}.${genome_ref}.vcf | sort -k1,1) | awk '{print $2"\t"$1}' > ${gagp_out_path}/${gagp_file_header}_discordant_ids_${genome_ref}_map.txt

	wc -l ${gagp_out_path}/${gagp_file_header}_discordant_ids_${genome_ref}_map.txt
	#overlap with original VCF
	join -eNA -o auto -t $'\t' -a 1 -11 <(sort -k1,1 ${gagp_out_path}/${gagp_file_header}_discordant_ids_${genome_ref}_map.txt) -21 <( awk '!($0 ~ /\#/) {print $1"_"$2"_"$4"_"$5"\t"$7}' ${gagp_out_path}/${gagp_file_header}.insertions.vcf | sort -k1,1)
	
done<${blacklist_work_dir}/gagp/gagp_meta_file.txt


#check that all ids are different are due to a filter being off
while read -a array
do
	gagp_file_header="${array[0]}"
	echo "${gagp_file_header}"
	
	#get panTro4 coordinates
	genome_ref=panTro4
	join -eNA -o auto -t $'\t' -a 1 -11 <(awk '{print}' ${blacklist_work_dir}/discord_files/novaseq_12_15_18_${gagp_file_header}_discordant_ids_${genome_ref}.txt | sort -k1,1) -21 <(awk '!($0 ~ /\#/) {print $3"\t"$1"_"$2"_"$4"_"$5}' ${hcondel_vcf_file} | sort -k1,1) | awk '{print $2"\t"$1}' > ${gagp_out_path}/novaseq_12_15_18_${gagp_file_header}_discordant_ids_${genome_ref}_map.temp.txt
	
	#overlap with unfiltered VCF coordinates
	join -eNA -o auto -t $'\t' -a 1 -11 <(sort -k1,1 ${gagp_out_path}/novaseq_12_15_18_${gagp_file_header}_discordant_ids_${genome_ref}_map.txt) -21 <( awk '!($0 ~ /\#/) {print $1"_"$2"_"$4"_"$5"\t"$7}' ${gagp_out_path}/${gagp_file_header}.${genome_ref}.normalized.vcf | sort -k1,1)
	
done<${blacklist_work_dir}/gagp/gagp_meta_file.txt


####### ##########

awk '{print $1}' /cluster_path_temp/hcondel_blacklist/Pan_paniscus_FINAL_INDEL.hcondel_overlapped.panTro4.txt | sort -k1,1 | uniq -c | awk '$1>1{print}'

grep 98005731 ${gagp_out_path}/Pan_troglodytes_FINAL_INDEL.hcondel_overlapped.panTro4.txt

grep 98005731 ${gagp_out_path}/Pan_troglodytes_FINAL_INDEL.hcondel_overlapped.panTro4.temp.txt

grep -E "97816236|97816237" ${gagp_out_path}/Pan_troglodytes_FINAL_INDEL.panTro4.normalized.vcf
grep -E "97816236|97816237" ${gagp_out_path}/Pan_troglodytes_FINAL_INDEL.hg38.normalized.vcf

grep -E "97816236|97816237" ${gagp_out_path}/Pan_troglodytes_FINAL_INDEL.panTro4.vcf
