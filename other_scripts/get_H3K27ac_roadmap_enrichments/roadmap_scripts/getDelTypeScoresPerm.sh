#!/bin/sh
#$ -cwd 
#$ -l h_vmem=10G
#$ -q long
#$ -V
#$ -j y
source /software_path/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/software_path/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/software_path/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/software_path/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/software_path/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/software_path/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib

permIter=$1
roadmapID=$2
modType=$3
cutOffVal=$4
posFile=$5
wigDir=$6
outDir=$7
permFileHeader=$8

#don't use chrSize for now- maybe use again in the future?


scriptdir=/cluster_path/ape_project/deletions_project/roadmap/scripts
#first gunzip the file
mkdir -p ${outDir}/${modType}

stat_out_file=${outDir}/${permFileHeader}_stat_out.txt
#echo "bigWigToBedGraph ${wiggleFileName} ${outDir}/${modType}/${roadmapID}_${modType}_${chr}_bedgraph.txt -chrom=${chr}"

#bigWigToBedGraph ${wiggleFileName} ${outDir}/${modType}/${roadmapID}_${modType}_${chr}_bedgraph.txt -chrom=${chr}
#I can choose to lift hg19 to hg38 here

#generate null sequences
#python /home/unix/sreilly/Scratch/TOOLS/kmersvm/scripts/nullseq_generate.py 
#intersect to get regions that lie in conserved area
#echo "bedtools intersect -a ${posFile} -b ${outDir}/${modType}/${roadmapID}_${modType}_${chr}_bedgraph.txt -loj > ${outDir}/${modType}/${roadmapID}_${modType}_chr_${chr}_intersected.txt"
bedtools intersect -sorted -a ${posFile} -b ${wigDir}/${roadmapID}-${modType}.pval.signal.bedgraph -loj > ${outDir}/${permFileHeader}_intersected.pval.txt
bedtools intersect -sorted -a ${posFile} -b ${wigDir}/${roadmapID}-${modType}.fc.signal.bedgraph -loj > ${outDir}/${permFileHeader}_intersected.fc.txt


#extract the scores after running bedtools interesect
#output will contain coordinate as well as average score/ if any part of region lies in a "peak" region
#echo "${scriptdir}/extractTypeScores_context.py ${outDir}/${modType}/${roadmapID}_${modType}_chr_${chr}_intersected.txt ${cutOffVal} ${outFile}"
python ${scriptdir}/extractTypeScores_context.py ${outDir}/${permFileHeader}_intersected.pval.txt ${cutOffVal} ${outDir}/${permFileHeader}_score_out.pval.txt
python ${scriptdir}/extractTypeScores_context.py ${outDir}/${permFileHeader}_intersected.fc.txt ${cutOffVal} ${outDir}/${permFileHeader}_score_out.fc.txt

#echo "python ${scriptdir}/extractTypeScores_context.py ${outDir}/${permFileHeader}_intersected.txt ${cutOffVal} ${score_out_file}"
#rm -f ${score_out_file}
#process the file and get stats
numSigCoordPval=$(awk '( ($(NF-1)==1)&& ($NF!=-1) ){print}' ${outDir}/${permFileHeader}_score_out.pval.txt | wc -l)
totalNumCoordPval=$(awk '$NF!=-1{print}' ${outDir}/${permFileHeader}_score_out.pval.txt | wc -l)
#propSigCoord=$(echo ${numSigCoord}/${totalNumCoord} | bc -l)
propSigCoordPval=$(echo -e "${numSigCoordPval}\t${totalNumCoordPval}" | awk -F"\t" '{print $1/$2}')
#echo ${score_out_file}
totalScorePval=$(cat ${outDir}/${permFileHeader}_score_out.pval.txt | awk '$4!=-1{ sum+=$4} END {print sum}' )
#echo ${totalScore}
totalLenPval=$(cat ${outDir}/${permFileHeader}_score_out.pval.txt | awk '$5!=-1{ sum+=$5} END {print sum}' )
#echo ${totalLen}
avgScorePval=$(echo -e "${totalScorePval}\t${totalLenPval}" | awk -F"\t" '{print $1/$2}')
#echo ${avgScore}



numSigCoordFc=$(awk '( ($(NF-1)==1)&& ($NF!=-1) ){print}' ${outDir}/${permFileHeader}_score_out.fc.txt | wc -l)
totalNumCoordFc=$(awk '$NF!=-1{print}' ${outDir}/${permFileHeader}_score_out.fc.txt | wc -l)
#propSigCoord=$(echo ${numSigCoord}/${totalNumCoord} | bc -l)
propSigCoordFc=$(echo -e "${numSigCoordFc}\t${totalNumCoordFc}" | awk -F"\t" '{print $1/$2}')
#echo ${score_out_file}
totalScoreFc=$(cat ${outDir}/${permFileHeader}_score_out.fc.txt | awk '$4!=-1{ sum+=$4} END {print sum}' )
#echo ${totalScore} 
totalLenFc=$(cat ${outDir}/${permFileHeader}_score_out.fc.txt | awk '$5!=-1{ sum+=$5} END {print sum}' )
#echo ${totalLen}
avgScoreFc=$(echo -e "${totalScoreFc}\t${totalLenFc}" | awk -F"\t" '{print $1/$2}')
#echo ${avgScore}


echo -e "${avgScorePval}\t${propSigCoordPval}\t${totalNumCoordPval}\t${avgScoreFc}\t${propSigCoordFc}\t${totalNumCoordFc}\t${roadmapID}"  > ${stat_out_file}

rm -f ${outDir}/${permFileHeader}_intersected.pval.txt
rm -f ${outDir}/${permFileHeader}_intersected.fc.txt
rm -f ${outDir}/${permFileHeader}_score_out.pval.txt
rm -f ${outDir}/${permFileHeader}_score_out.fc.txt
