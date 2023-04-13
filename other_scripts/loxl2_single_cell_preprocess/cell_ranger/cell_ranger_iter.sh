id=$1
transcriptome=$2
fastqs=$3
sample=$4
expect_cells=$5
localcores=$6
localmem=$7


echo "cellranger count --id=${id} --transcriptome=${transcriptome} --fastqs=${fastqs} --sample=${sample} --expect-cells=${expect_cells} --localcores=${localcores} --localmem=${localmem}"
cellranger count --id=${id} \
                   --transcriptome=${transcriptome}\
                   --fastqs=${fastqs} \
                   --sample=${sample} \
                   --expect-cells=${expect_cells} \
                   --localcores=${localcores} \
                   --localmem=${localmem}
