config=$1
fastqR1=$2
fastqR2=$3
whitelist=$4
outdir=$5
log=$6
sample=$7

IronThrone-GoT --run linear --config ${config} --thread 1 --umilen 12 --verbose --fastqR1 ${fastqR1} --fastqR2 ${fastqR2} --whitelist ${whitelist} --outdir ${outdir} --log ${log} --sample ${sample}

