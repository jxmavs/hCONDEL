
#obtain file list names
#ls -l /cluster_path_temp/del_encode/coord_split/${real_data_out_header}_all_coord_extended_with_fimo_id_sorted.fa.* > ${out_path}/all_fimo_file_names.txt
#combine all files from fimo run into one master file
out_header=$1
#the file header for the resulting output
file_out_header=$2
out_path=$3
pairwise_comp_file=$4

rm -f ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.txt
rm -f ${out_path}/${out_header}_fimo_cleaned_sum_diff_output_final.txt


######################## for max diff ################################# 
#add header

preHeaderInfo="${file_out_header}_max_diff"
headerLine="SeqName\t${preHeaderInfo}_TF"
#check that the order is correct
species_arr=("panTro4" "hg38" "rheMac8")
meta_data_col_names=("score" "startPos" "endPos" "strand" "pval" "overlapDel")

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

echo -e "${headerLine}" > ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.txt


 
#################### for sum diff ###################################

preHeaderInfo="${file_out_header}_sum_diff"
headerLine="SeqName\t${preHeaderInfo}_TF"
species_arr=("panTro4" "hg38" "rheMac8")
meta_data_col_names=("avg_score")

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
 
headerLine="${headerLine}\t${preHeaderInfo}_TF_annote"
echo -e "${headerLine}" > ${out_path}/${out_header}_fimo_cleaned_sum_diff_output_final.txt

############################

while read -a array
do
    
    input_file=${array[0]}
    suffix=$( echo ${input_file} | awk -F"." '{print $NF}' )
    cat ${out_path}/${out_header}_fimo_cleaned_max_diff_output.${suffix}.txt >>  ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.txt

    cat ${out_path}/${out_header}_fimo_cleaned_sum_diff_output.${suffix}.txt >>  ${out_path}/${out_header}_fimo_cleaned_sum_diff_output_final.txt

done<${out_path}/all_fimo_file_names.txt


######## join all files together ################ 

join -11 ${out_path}/${out_header}_fimo_cleaned_max_diff_output_final.txt -21 ${out_path}/${out_header}_fimo_cleaned_sum_diff_output_final.txt | awk 'BEGIN{OFS="\t";} {$1=$1}1' > ${out_path}/${out_header}_fimo_cleaned_output_final.txt 



