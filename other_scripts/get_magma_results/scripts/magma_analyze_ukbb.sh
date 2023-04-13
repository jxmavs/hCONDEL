#similar to magma_analyze.sh, except with ukbb data

phenotype_code=$1
download_link=$2
download_file_header=$3
out_path=$4
del_set_header=$5
geneloc_file=$6
refdata=$7
set_file=$8
additional_meta_file=$9

workdir=${out_path}/${phenotype_code}
rm -fr ${workdir}
mkdir -p ${workdir}

#download file

wget -nv ${download_link} -O ${workdir}/${download_file_header}.bgz
mv ${workdir}/${download_file_header}.bgz ${workdir}/${download_file_header}.gz

#process file

analysis_header=${phenotype_code}_magma_analysis
snploc_file_initial=${workdir}/${download_file_header}.gz

#generate snploc file
gunzip -c ${snploc_file_initial} | awk 'NR>1{print $1}' | awk -F":" '{print $1"_"$2"_"$3"_"$4"\t"$1"\t"$2}' > ${workdir}/${analysis_header}_magma_snploc_input.txt
#generate pval_file
gunzip -c ${snploc_file_initial} | awk 'NR>1{print $NF"\t"$(NF-6)}' > ${workdir}/${analysis_header}_magma_pval_input.2.txt
paste <(awk '{print $1}' ${workdir}/${analysis_header}_magma_snploc_input.txt) <(cat ${workdir}/${analysis_header}_magma_pval_input.2.txt) > ${workdir}/${analysis_header}_magma_pval_input.txt
sed -i 's#NaN#NA#g' ${workdir}/${analysis_header}_magma_pval_input.txt
rm -f ${workdir}/${analysis_header}_magma_pval_input.2.txt
rm -f ${workdir}/${download_file_header}.gz

#a header is not allowed, file must have three columns containing the SNP ID, chromosome and base pair position
snploc_file=${workdir}/${analysis_header}_magma_snploc_input.txt

annot_prefix=${workdir}/${analysis_header}_${del_set_header}_annot

#file with previously computed SNP p-values (in columns 'SNP' and 'P' in the file [PVAL_FILE])
pval_file=${workdir}/${analysis_header}_magma_pval_input.txt
gene_prefix=${workdir}/${analysis_header}_${del_set_header}_gene

gs_prefix=${workdir}/${analysis_header}_${del_set_header}_gs

#get N, should only be one N
awk '{print $3}' ${workdir}/${analysis_header}_magma_pval_input.txt | sort | uniq > ${workdir}/${analysis_header}_${del_set_header}_N.txt  

#run magma

magma --annotate --snp-loc ${snploc_file} --gene-loc ${geneloc_file} --out ${annot_prefix}

magma --bfile ${refdata} --pval ${pval_file} ncol=3 --gene-annot ${annot_prefix}.genes.annot --out ${gene_prefix}

magma --gene-results ${gene_prefix}.genes.raw --set-annot ${set_file} col=1,2 --gene-covar ${additional_meta_file} filter-read=cons_num,cons_prop --model condition=cons_num,cons_prop --out ${gs_prefix}_condition_cons_prop

#remove files
rm -f ${workdir}/${analysis_header}_magma_snploc_input.txt
rm -f ${workdir}/${analysis_header}_magma_pval_input.txt
