#combine all files from fimo run into one master file
out_header=$1
#the file header for the resulting output
file_out_header=$2
out_path=$3
pairwise_comp_file=$4
tf_name=$5
seqFileNamesFile=$6


######################## for max diff ################################# 
#add header

preHeaderInfo="${file_out_header}_max_diff"
headerLine="SeqName\t${preHeaderInfo}_TF"
#check that the order is correct
species_arr=("panTro4" "hg38" "rheMac8")
meta_data_col_names=("score" "startPos" "endPos" "strand" "pval" "overlapDel")


rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.${tf_name}.txt

for species in "${species_arr[@]}"
do
	for meta_col_name in "${meta_data_col_names[@]}"
	do
		headerIter="${preHeaderInfo}_${species}_${meta_col_name}"
		headerLine="${headerLine}\t${headerIter}"
	done
	
done

while read -a array
do
	species1="${array[0]}"
	species2="${array[1]}"
	headerIter="${preHeaderInfo}_${species2}_${species1}_diff"
	headerLine="${headerLine}\t${headerIter}"
done<${pairwise_comp_file}

echo -e "${headerLine}" > ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.${tf_name}.txt

############################

while read -a array
do
	
	input_file=${array[0]}
	file_id_suffix=$( echo ${input_file} | awk -F"." '{print $NF}' )
	suffix=${file_id_suffix}"."${tf_name}
	
	cat ${out_path}/${out_header}_fimo_cleaned_max_diff_output.${suffix}.txt >> ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.${tf_name}.txt

done<${seqFileNamesFile}



