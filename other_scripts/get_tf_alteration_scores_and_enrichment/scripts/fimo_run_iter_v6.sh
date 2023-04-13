#v5 takes encodeIDMapFileName and maps file ID to TF name 
#10/23/21 corrected typo with not properly acccounting for hg38 deletion position, changed from ${scriptdir}/fimo_get_max_diff_v2.py to ${scriptdir}/fimo_get_max_diff_v3.py, changed from ${scriptdir}/fimo_get_sum_diff_v2.py to ${scriptdir}/fimo_get_sum_diff_v3.py
input_file=$1
out_header=$2
motif_file=$3
pairwise_comp_file=$4
lookAtCenteredOnly=$5
cutOffVal=$6
useEncodeAnnotation=$7
encodeFileName=$8
encodeIDMapFileName=$9
lookAtExpressedTFs=${10}
expressedTFsFileName=${11}
memeFileName=${12}
out_path=${13}

scriptdir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts
python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin

#get the file suffix
suffix=$( echo ${input_file} | awk -F"." '{print $NF}' )

#run fimo
fimo --thresh 1 --verbosity 1 --no-qvalue --text ${motif_file} ${input_file} > ${out_path}/${out_header}_fimo_output.${suffix}.txt 

#map hocomoco ID to genecard ID
analysis_path=/cluster_path/ape_project/deletions_project/meme_suite_analysis

join -t $'\t' -11 -21 <(sort -k1,1 ${out_path}/${out_header}_fimo_output.${suffix}.txt)  <(sort -k1,1 ${memeFileName} | awk '{print $1"\t"$2}') > ${out_path}/${out_header}_fimo_output.${suffix}.temp.txt

paste <( awk -F"\t" '{print $NF}' ${out_path}/${out_header}_fimo_output.${suffix}.temp.txt )  <(cut -d$'\t' -f3-8,10 ${out_path}/${out_header}_fimo_output.${suffix}.temp.txt ) >  ${out_path}/${out_header}_fimo_output.${suffix}.withTFId.txt

mv ${out_path}/${out_header}_fimo_output.${suffix}.withTFId.txt ${out_path}/${out_header}_fimo_output.${suffix}.txt

#after running FIMO, get the max score difference between the max of species 2 and species 1
${python_dir}/python ${scriptdir}/fimo_get_max_diff_v3.py ${out_path}/${out_header}_fimo_output.${suffix}.txt ${pairwise_comp_file} ${out_path}/${out_header}_fimo_max_diff_output.${suffix}.txt ${lookAtCenteredOnly} 


#after getting max difference, check to see if either human or chimp has pval<cutoff, and if a coordinate has multiple TFs, keep the one with the max difference
#mark the TFs with annotation also if desired 
#if useEncodeAnnotation is set, then only look at the ones with encode annotation

#the species p value indeces, change if we change the output of the maxdiff output from fimo_get_max_diff.py
#species1 is human, species 2 is chimp
#we ignore macaque for now - the same with looking to see if TF is expressed - just check human and chimp for now

species1pvalind=6
species2pvalind=12
maxDiffInd=20

${python_dir}/python ${scriptdir}/parse_max_diff_v4.py ${out_path}/${out_header}_fimo_max_diff_output.${suffix}.txt ${species1pvalind} ${species2pvalind} ${maxDiffInd} ${cutOffVal} ${useEncodeAnnotation} ${encodeFileName} ${encodeIDMapFileName} ${lookAtExpressedTFs} ${expressedTFsFileName} ${memeFileName} ${out_path}/${out_header}_fimo_cleaned_max_diff_output.${suffix}.txt

########## get sum difference ###################

${python_dir}/python ${scriptdir}/fimo_get_sum_diff_v3.py ${out_path}/${out_header}_fimo_output.${suffix}.txt ${pairwise_comp_file} ${out_path}/${out_header}_fimo_sum_diff_output.${suffix}.txt ${lookAtCenteredOnly} ${cutOffVal}

sumDiffInd=5

${python_dir}/python ${scriptdir}/parse_sum_diff_v3.py ${out_path}/${out_header}_fimo_sum_diff_output.${suffix}.txt ${sumDiffInd} ${useEncodeAnnotation} ${encodeFileName} ${encodeIDMapFileName} ${lookAtExpressedTFs} ${expressedTFsFileName} ${memeFileName} ${out_path}/${out_header}_fimo_cleaned_sum_diff_output.${suffix}.txt



rm -f ${out_path}/${out_header}_fimo_output.${suffix}.txt
rm -f ${out_path}/${out_header}_fimo_max_diff_output.${suffix}.txt
rm -f ${out_path}/${out_header}_fimo_sum_diff_output.${suffix}.txt 

rm -f ${out_path}/${out_header}_fimo_output.${suffix}.temp.txt




 
