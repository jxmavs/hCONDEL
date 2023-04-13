#This script generates permutations for ${modType} across all cell types and checks for enrichment

#analyze it one chromosome at a time
roadmap_dir=/cluster_path/ape_project/deletions_project/roadmap
wrkdir=${roadmap_dir}/work

modType=$1
roadmap_ID_file=$2

#cutOffVal will be the threshold value of what we consider significant or not
cutOffVal=$3
outDir=$4
wigDir=$5
out_header=$6
posFile=$7
#seperate out deletions into files based on individual chromosomes for speedup

#modType=H3K27ac
#type_dir=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${modType}
#roadmap_ID_file=/cluster_path/ape_project/deletions_project/roadmap/nonimputed/${modType}/${modType}_file_list.txt
#cutOffVal=2
#out_header=del_roadmap_nonimputed_${modType}

ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human

#first iterate over cell type, 
#then iterate over chrmosome
rm -f ${wrkdir}/${out_header}_${modType}_all_runs.txt 
while IFS=$'\t' read -a array
do
    roadmapID=${array[0]}
    #while read -a array2
    #do
    #chr=${array2[0]}
    #chrSize=$(awk -v chrID=$chr '$1==chrID {print $2}' ${ref_dir}/hg19.chrom.sizes)
    #$posFile are the null sequences generated from gen_null_seq_submit.sh
    outFileHeader=${out_header}_${roadmapID}_${modType}_cutOff_${cutOffVal}
    echo "bash ${roadmap_dir}/scripts/getDelTypeScoresRealData.sh ${roadmapID} ${modType} ${cutOffVal} ${posFile} ${wigDir} ${outDir} ${outFileHeader}" >> ${wrkdir}/${out_header}_${modType}_all_runs.txt 
    #done<${wrkdir}/chr_list.txt

done<${roadmap_ID_file}
echo -e "${wrkdir}/${out_header}_${modType}_all_runs.txt"
#run everything

del_dir=/cluster_path/ape_project/deletions_project
totalMem=1G
totalNumTasks=$(wc -l ${wrkdir}/${out_header}_${modType}_all_runs.txt  | awk '{print $1}')   
#sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${wrkdir}/${out_header}_${modType}_getTypeScoresAll.sh
totalTime=30:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${wrkdir}/${out_header}_${modType}_getTypeScoresAll.sh
rm -f ${wrkdir}/${out_header}_${modType}.log 
qsub -o ${wrkdir}/${out_header}_${modType}.log ${wrkdir}/${out_header}_${modType}_getTypeScoresAll.sh ${wrkdir}/${out_header}_${modType}_all_runs.txt

#then merge files afterwards using merge_roadmap_score_files.sh
#final output is bed position, ID, then indivscores in region, then avg score in region
