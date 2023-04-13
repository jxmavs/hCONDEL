snploc_file=$1
geneloc_file=$2
annot_prefix=$3
refdata=$4
pval_file=$5
gene_prefix=$6
N=$7
set_file=$8
additional_meta_file=$9
gs_prefix=${10}

magma --annotate --snp-loc ${snploc_file} --gene-loc ${geneloc_file} --out ${annot_prefix}

#if N is not numeric, then use file 
re='^[0-9]+$'
if [[ $N =~ $re ]] ; then
	magma --bfile ${refdata} --pval ${pval_file} N=${N} --gene-annot ${annot_prefix}.genes.annot --out ${gene_prefix}
else
	magma --bfile ${refdata} --pval ${pval_file} ncol=3 --gene-annot ${annot_prefix}.genes.annot --out ${gene_prefix}
fi

if [[ ${additional_meta_file} == "NULL" ]]; then
    magma --gene-results ${gene_prefix}.genes.raw --set-annot ${set_file} col=1,2 --out ${gs_prefix}
else
	magma --gene-results ${gene_prefix}.genes.raw --set-annot ${set_file} col=1,2 --gene-covar ${additional_meta_file} filter-read=cons_num,cons_prop --out ${gs_prefix}

	magma --gene-results ${gene_prefix}.genes.raw --set-annot ${set_file} col=1,2 --gene-covar ${additional_meta_file} filter-read=cons_num,cons_prop --model condition=cons_num,cons_prop --out ${gs_prefix}_condition_cons_prop
fi
