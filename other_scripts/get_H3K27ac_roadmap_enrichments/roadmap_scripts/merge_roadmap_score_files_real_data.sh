
#this script merges output score files from all roadmap cell lines

modType=$1
roadmap_ID_file=$2
cutOffVal=$3
out_header=$4
outDir=$5

statOutFile=${outDir}/all_${out_header}_${modType}_cutOff_${cutOffVal}_stat_out.txt
echo -e "${statOutFile}"
rm -f ${statOutFile}

while IFS=$'\t' read -a array
do
    #simply add on cell type name
    roadmapID=${array[0]}
    cellType=${array[1]}
    #echo ${cellType}
    iterOutFileName=${outDir}/${out_header}_${roadmapID}_${modType}_cutOff_${cutOffVal}_stat_out.txt
    paste <(awk '{print}' ${iterOutFileName}) <( echo ${cellType} ) >> ${statOutFile}

done<${roadmap_ID_file}


