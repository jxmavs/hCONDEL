#!/bin/sh
#v2 incoporates logOdds as well
permInd=$1
out_header=$2
out_path=$3
cons_seq_comp_stat_file=$4
orig_seq_comp_stat_file=$5
lengthBoundPct=$6
misMatchBoundPct=$7
GCBoundPct=$8
logOddsBoundPct=$9
sampleNullSeqMaxTrys=${10}
sampleNullSeqConsDevPct=${11}


perm_dir=/cluster_path_temp/roadmap/permutations

crossSpeciesScriptDir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin

#python /home/unix/sreilly/Scratch/TOOLS/kmersvm/scripts/nullseq_generate.py -x 1 -r ${permInd} -R -G -e ${exclude_coord} -o ${outdir}/${permInd}_${out_header}_list.txt ${input_coord} hg19 ${perm_dir}

#sort -k1,1 -k2,2n ${outdir}/${permInd}_${out_header}_list.txt | awk '{print $0"\t"$1"|"$2"|"$3}' > ${outdir}/${permInd}_${out_header}_list_sorted.txt
#mv ${outdir}/${permInd}_${out_header}_list_sorted.txt ${outdir}/${permInd}_${out_header}_list.txt

permFileHeader=${permInd}_${out_header}

numSeqOrigFile=$(wc -l ${orig_seq_comp_stat_file} | awk '{print $1}')
${python_dir}/python ${crossSpeciesScriptDir}/genRandSampleV4.py ${cons_seq_comp_stat_file} ${orig_seq_comp_stat_file} ${numSeqOrigFile} ${lengthBoundPct} ${misMatchBoundPct} ${GCBoundPct} ${logOddsBoundPct} ${out_path}/${permFileHeader}_cons_sequences_info.txt
#genRandSampleV2 outputs statistics (the amount of deviation from the original set) with each sampled id


#add one to human deletion position here
paste <( awk '{print $1}' ${out_path}/${permFileHeader}_cons_sequences_info.txt |  awk -F'[|]' '{print $4"\t"$5"\t"$6}' ) <( awk '{print $NF}'  ${orig_seq_comp_stat_file} | awk -F'[|]' '{print $4"\t"$5"\t"$6"\t"$10"\t"$11"\t"$12+1"\t"$0}' ) > ${out_path}/${permFileHeader}_cons_sequences_input.txt 

#input contains fields randomly sampled conserved position, original conserved position, original deletion position, id
${python_dir}/python ${crossSpeciesScriptDir}/sampleNullSeqRandConsBackgroundHuman.py ${out_path}/${permFileHeader}_cons_sequences_input.txt ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${out_path}/${permFileHeader}_list.txt

sort -k1,1 -k2,2n ${out_path}/${permFileHeader}_list.txt > ${out_path}/${permFileHeader}_list_sorted.txt
mv ${out_path}/${permFileHeader}_list_sorted.txt ${out_path}/${permFileHeader}_list.txt

#rm -f ${out_path}/${permFileHeader}_cons_sequences_input.txt
#rm -f ${out_path}/${permFileHeader}_cons_sequences_info.txt
#python ${crossSpeciesScriptDir}/genRandSampleV3.py ${cons_seq_comp_stat_file} ${orig_seq_comp_stat_file} ${numSeqOrigFile} ${lengthBoundPct} ${misMatchBoundPct} ${out_path}/${permFileHeader}_cons_sequences_info.txt 
#genRandSampleV2 outputs statistics (the amount of deviation from the original set) with each sampled id
#modify file for input into sampleNullSeq.py

#awk '{print $1}' ${out_path}/${permFileHeader}_cons_sequences_info.txt | awk -F'|' '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$0}' > ${out_path}/${permFileHeader}_cons_sequences_input.txt

#python ${crossSpeciesScriptDir}/sampleNullSeq.py ${out_path}/${permFileHeader}_cons_sequences_input.txt ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${out_path}/${permFileHeader}_list.txt

