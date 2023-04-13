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

permSequencesFinalFile=${out_path}/${permFileHeader}_list.final.txt
rm -f ${permSequencesFinalFile}

finishedSampling=0

#if conserved blocks are already drawn, no need to redraw
initialConsDrawn=0
#to keep track of previously drawn sequences
prevDrawnInfoFileName="NA"

ref_genome_chr_sizes_path=/cluster_path/ape_project/deletions_project/reference_genome_files/human/hg38.chrom.sizes

consSequenceInputFile=${out_path}/${permFileHeader}_cons_sequences_input.txt

#keep track of previously drawn sequences to not resample those in while loop
rm -f ${out_path}/${permFileHeader}.orig_input_info.prev_drawn.txt
#look at unmapped seq
rm -f ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.temp.txt

#remove additional files that are appended to just in case they previously existed
rm -f ${out_path}/${permFileHeader}_human_coord.temp.txt
rm -f ${out_path}/${permFileHeader}_chimp_coord.temp.txt
rm -f ${out_path}/${permFileHeader}_log.txt

while [ ${finishedSampling} -eq 0 ]
do
	numSeqOrigFile=$(wc -l ${orig_seq_comp_stat_file} | awk '{print $1}')
	
	if [ ${initialConsDrawn} -eq 0 ]
	then
		${python_dir}/python ${crossSpeciesScriptDir}/genRandSampleV4.py ${cons_seq_comp_stat_file} ${orig_seq_comp_stat_file} ${numSeqOrigFile} ${lengthBoundPct} ${misMatchBoundPct} ${GCBoundPct} ${logOddsBoundPct} ${out_path}/${permFileHeader}_cons_sequences_info.txt
	
		#genRandSampleV2 outputs statistics (the amount of deviation from the original set) with each sampled id

		#10/22/21 revised:put chimp deletion coordinates last, it was $10"\t"$11"\t"$12+1, which was used to get 1 bp for human seq sampling, now its $7"\t"$8"\t"$9
		#output is sampled human cons cord, original human cons cord, chimp del coord, deletion id
		paste <( awk '{print $1}' ${out_path}/${permFileHeader}_cons_sequences_info.txt |  awk -F'[|]' '{print $4"\t"$5"\t"$6}' ) <( awk '{print $NF}'  ${orig_seq_comp_stat_file} | awk -F'[|]' '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$0}' ) > ${out_path}/${permFileHeader}_cons_sequences_input.txt
	fi
	
	#input contains fields randomly sampled conserved position, original conserved position, original deletion position, id
	${python_dir}/python ${crossSpeciesScriptDir}/sampleNullSeqRandConsBackgroundHuman_past_cons_block_v2.py ${consSequenceInputFile} ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${out_path}/${permFileHeader}_list.txt ${ref_genome_chr_sizes_path} ${prevDrawnInfoFileName}
	
	cat ${out_path}/${permFileHeader}_list.txt | sort -k1,1 -k2,2n > ${out_path}/${permFileHeader}_list_sorted.txt
	mv ${out_path}/${permFileHeader}_list_sorted.txt ${out_path}/${permFileHeader}_list.txt

	#rm -f ${out_path}/${permFileHeader}_cons_sequences_input.txt
	#rm -f ${out_path}/${permFileHeader}_cons_sequences_info.txt
	#python ${crossSpeciesScriptDir}/genRandSampleV3.py ${cons_seq_comp_stat_file} ${orig_seq_comp_stat_file} ${numSeqOrigFile} ${lengthBoundPct} ${misMatchBoundPct} ${out_path}/${permFileHeader}_cons_sequences_info.txt 
	#genRandSampleV2 outputs statistics (the amount of deviation from the original set) with each sampled id
	#modify file for input into sampleNullSeq.py

	#awk '{print $1}' ${out_path}/${permFileHeader}_cons_sequences_info.txt | awk -F'|' '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$0}' > ${out_path}/${permFileHeader}_cons_sequences_input.txt

	#python ${crossSpeciesScriptDir}/sampleNullSeq.py ${out_path}/${permFileHeader}_cons_sequences_input.txt ${sampleNullSeqMaxTrys} ${sampleNullSeqConsDevPct} ${out_path}/${permFileHeader}_list.txt
	
	awk '{print $1"\t"$2"\t"$3"\t"$4"#"$1"|"$2"|"$3"\t0\t+"}' ${out_path}/${permFileHeader}_list.txt > ${out_path}/${permFileHeader}_list.revised.txt
	#lift to panTro4, if deleted, or doesn't match exact sequence, then sample again
	liftChainFile=/cluster_path/ape_project/deletions_project/liftover_files/hg38ToPanTro4.over.chain.gz
	liftOver -minMatch=1 ${out_path}/${permFileHeader}_list.revised.txt ${liftChainFile} ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.txt ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt

	#human coord first/chimp second, followed by id at end
	paste <(awk '{print $4}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.txt | awk -F"#" '{print $2}' | awk -F"|" '{print $1"\t"$2"\t"$3}') <(awk '{print $1"\t"$2"\t"$3}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.txt) <(awk '{print $4}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.txt | awk -F"#" '{print $1}') <(awk '{print $6}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.txt) > ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.txt


	#check to see if any coordinates overlap, if it does, put it aside to discard, check for both human and chimp
	#add on already drawn coordinates
	if [[ -s ${permSequencesFinalFile} ]]
	then
		awk '{print $1"\t"$2"\t"$3"\t"$7}' ${permSequencesFinalFile} > ${out_path}/${permFileHeader}_human_coord.temp.txt
		awk '{print $4"\t"$5"\t"$6"\t"$7}' ${permSequencesFinalFile} > ${out_path}/${permFileHeader}_chimp_coord.temp.txt
	fi
	
	#add on coordinates in this iteration
	awk '{print $1"\t"$2"\t"$3"\t"$7}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.txt >> ${out_path}/${permFileHeader}_human_coord.temp.txt
	awk '{print $4"\t"$5"\t"$6"\t"$7}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.txt >> ${out_path}/${permFileHeader}_chimp_coord.temp.txt
	
	sort -k1,1 -k2,2n ${out_path}/${permFileHeader}_human_coord.temp.txt > ${out_path}/${permFileHeader}_human_coord.temp.sorted.txt
	mv ${out_path}/${permFileHeader}_human_coord.temp.sorted.txt ${out_path}/${permFileHeader}_human_coord.temp.txt

	sort -k1,1 -k2,2n ${out_path}/${permFileHeader}_chimp_coord.temp.txt > ${out_path}/${permFileHeader}_chimp_coord.temp.sorted.txt
	mv ${out_path}/${permFileHeader}_chimp_coord.temp.sorted.txt ${out_path}/${permFileHeader}_chimp_coord.temp.txt
		
	bedtools intersect -sorted -loj -a ${out_path}/${permFileHeader}_human_coord.temp.txt -b ${out_path}/${permFileHeader}_human_coord.temp.txt | awk '{print $4}' | uniq -c | awk '$1>1{print $2}' > ${out_path}/${permFileHeader}_human_coord.temp.intersect.ids
	bedtools intersect -sorted -loj -a ${out_path}/${permFileHeader}_chimp_coord.temp.txt -b ${out_path}/${permFileHeader}_chimp_coord.temp.txt | awk '{print $4}' | uniq -c | awk '$1>1{print $2}' > ${out_path}/${permFileHeader}_chimp_coord.temp.intersect.ids

	#combine human + chimp intersected ids
	cat ${out_path}/${permFileHeader}_human_coord.temp.intersect.ids > ${out_path}/${permFileHeader}_all.temp.intersect.ids
	cat ${out_path}/${permFileHeader}_chimp_coord.temp.intersect.ids >> ${out_path}/${permFileHeader}_all.temp.intersect.ids	

	#get the ids from this iteration of the loop
	sort -k1,1 ${out_path}/${permFileHeader}_all.temp.intersect.ids | uniq > ${out_path}/${permFileHeader}_all.temp.intersect.uniq.ids 

	#get cleaned file
	join -v2 -11 <(sort -k1,1 ${out_path}/${permFileHeader}_all.temp.intersect.uniq.ids) -27 <(sort -k7,7 ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.txt) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$8}' > ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.intersect_removed.txt	
	
	if [[ ! -s ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt ]] && [[ ! -s ${out_path}/${permFileHeader}_all.temp.intersect.uniq.ids ]] 
	then		
		finishedSampling=1
		#human coord first/chimp second, add onto final file
		cat ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.intersect_removed.txt >> ${permSequencesFinalFile}
	else
		initialConsDrawn=1
		#keep track of unmapped for debugging purposes
		cat ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt >> ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.temp.txt

        cat ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.intersect_removed.txt >> ${permSequencesFinalFile}

		#use the same conserved sequence sample another deletion
		awk 'NR%2==0{print}' ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt | awk '{print $4}' | awk -F"#" '{print $1}' | sort -k1,1 > ${out_path}/${permFileHeader}_unmapped.ids.txt
		#join -11 ${out_path}/${permFileHeader}_unmapped.ids.txt -210 <( sort -k10,10 ${out_path}/${permFileHeader}_cons_sequences_input.txt ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$1}' > ${out_path}/${permFileHeader}.cons_stats.txt
	
		wc -l ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt  >> ${out_path}/${permFileHeader}_log.txt	
		wc -l ${consSequenceInputFile} >> ${out_path}/${permFileHeader}_log.txt
		
		#for debugging, copy and compare later, rearrange
		cp ${consSequenceInputFile} ${out_path}/${permFileHeader}.gen_rand_empirical_data_input.txt
		
		#for resampling
		cat ${out_path}/${permFileHeader}_unmapped.ids.txt > ${out_path}/${permFileHeader}_resample.ids.txt
		cat ${out_path}/${permFileHeader}_all.temp.intersect.uniq.ids >> ${out_path}/${permFileHeader}_resample.ids.txt
		
		join -11 <(sort -k1,1 ${out_path}/${permFileHeader}_resample.ids.txt) -210 <(sort -k10,10 ${out_path}/${permFileHeader}.gen_rand_empirical_data_input.txt) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$1}' > ${out_path}/${permFileHeader}.orig_input_info.txt
		
		#add onto already previously drawn sequences 
		join -11 <(awk '{print $7"\t"$1"\t"$2"\t"$3}' ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.intersect_removed.txt | sort -k1,1) -210 <(sort -k10,10 ${out_path}/${permFileHeader}.gen_rand_empirical_data_input.txt) | awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$2"\t"$3"\t"$4"\t"$1}' >> ${out_path}/${permFileHeader}.orig_input_info.prev_drawn.txt	
		
		rm -f ${out_path}/${permFileHeader}_list.revised.panTro4.unmapped.txt 
		
		consSequenceInputFile=${out_path}/${permFileHeader}.orig_input_info.txt
		prevDrawnInfoFileName=${out_path}/${permFileHeader}.orig_input_info.prev_drawn.txt
	
		wc -l ${out_path}/${permFileHeader}_human_coord.temp.intersect.ids >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${out_path}/${permFileHeader}_chimp_coord.temp.intersect.ids >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${out_path}/${permFileHeader}_all.temp.intersect.uniq.ids >> ${out_path}/${permFileHeader}_log.txt	
		wc -l ${out_path}/${permFileHeader}_list.revised.panTro4.mapped.temp.intersect_removed.txt >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${permSequencesFinalFile} >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${out_path}/${permFileHeader}_unmapped.ids.txt >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${out_path}/${permFileHeader}_resample.ids.txt >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${consSequenceInputFile} >> ${out_path}/${permFileHeader}_log.txt
		wc -l ${prevDrawnInfoFileName} >> ${out_path}/${permFileHeader}_log.txt
	fi	

done

sort -k7,7 ${permSequencesFinalFile} > ${permSequencesFinalFile}.sorted
mv ${permSequencesFinalFile}.sorted ${permSequencesFinalFile}
 
#below script is adapted from cell_rna_process_commands.sh for FIMO sampling
awk '{print $1"\t"$2"\t"$3"\t+\t"$7}' ${permSequencesFinalFile} > ${out_path}/${permFileHeader}_human_coord.txt
awk '{print $4"\t"$5"\t"$6"\t"$8"\t"$7}' ${permSequencesFinalFile} > ${out_path}/${permFileHeader}_chimp_coord.txt

#set to 60 because largest tf motif is 25ish bp in JASPAR
maxBPWindow=60
tf_enrichment_files_dir=/cluster_path/ape_project/deletions_project/tf_enrichment/files
species_ref_meta_file=${tf_enrichment_files_dir}/species_ref_input_file.txt
scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts

#rm -f ${out_path}/${permFileHeader}_all_coord_extended_with_fimo_id.fa

rm -f ${out_path}/${permFileHeader}_all_fastas_list.txt
while read -a array
do
    speciesName="${array[0]}"
    ref_genome="${array[1]}"
	chrSizeFile="${array[2]}"
	
	#extend coordinates first
	${python_dir}/python ${scriptdir}/extendCoordinatesCHIPV3b.py ${out_path}/${permFileHeader}_${speciesName}_coord.txt ${chrSizeFile} ${maxBPWindow} ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt
	
	#don't sort because sequences can be scrambled in matching		
	#sort -k1,1 -k2,2n ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt > ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_temp.txt
	#mv ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_temp.txt ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt

	#get fasta
	fasta_ref_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/${speciesName}

	#create id, add on the species, along with the extended coordinates and the centered coordinate prior to extension
	paste -d"|" <( awk -v ref_genome=${ref_genome} '{print $7"#"ref_genome}' ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt ) <( cut -f1-6 ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt | awk 'BEGIN{OFS="|";} {$1=$1}1'  ) > ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_fimo_id.txt

	#revise file header for coordinate
	paste <(cut -f1-3 ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_bp_b.txt ) <(cat ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_fimo_id.txt)  > ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.txt
	bedtools getfasta -name -fi ${fasta_ref_dir}/${ref_genome}.fa -bed ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.txt -fo ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa

	if [[ "${speciesName}" == "human" ]]
	then
		#for human, remove deleted portion in sequence, revise coordinate to account for removal
		python ${scriptdir}/removeDelHumanSeq.py ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.del_seq_removed.fa
		#keep track for debugging
		cp ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.pre_del_seq_removed.fa
		mv ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.del_seq_removed.fa ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa
	fi
	
	#cat ${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa >> ${out_path}/${permFileHeader}_all_coord_extended_with_fimo_id.fa 
	echo "${out_path}/${permFileHeader}_${speciesName}_coord_extended_${maxBPWindow}_with_fimo_id.fa" >> ${out_path}/${permFileHeader}_all_fastas_list.txt
done<${species_ref_meta_file}

#sort, don't use because possible that duplicate hcondels can be sampled, see maxTrys
#awk '!/^>/ { next } { getline seq } { print $0"\t"seq }'${out_path}/${permFileHeader}_all_coord_extended_with_fimo_id.fa | sort -k1,1 | awk '{print $1"\n"$2}' > ${out_path}/${permFileHeader}_all_coord_extended_with_fimo_id_sorted.fa

#don't use, because coord id is matching the original hcondel, but can: combine seperately, as worried about duplicate sequences being sampled if I simply combine, then sort 
python ${scriptdir}/combineSpeciesFastas.py ${out_path}/${permFileHeader}_all_fastas_list.txt ${out_path}/${permFileHeader}_all_coord_extended_with_fimo_id.fa
