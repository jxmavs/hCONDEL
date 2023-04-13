#aligns human to chimp sequence to get mismatch stats
out_path_mismatch_analysis=$1
header_name=$2
seq_name=$3
chimp_seq=$4
human_seq=$5
seq_id=$6

python_dir=/software_path/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin
script_dir=/cluster_path/ape_project/deletions_project/cross_species_chip/scripts
#create fasta from table 
echo -e "${seq_name}\n${chimp_seq}" > ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa

echo -e "${seq_name}\n${human_seq}" > ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.fa

#blast
#create database
#added /dev/null 4/23/19

makeblastdb -in ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa -dbtype nucl -logfile "/dev/null"

#dust is the low complexity filter
blastn -db ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa -query ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.fa -out ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats -penalty -3 -reward 2 -gapopen 5 -gapextend 2 -dust no -word_size 10 -evalue 1 -outfmt "6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch gapopen gaps evalue bitscore"

#sort by bit score (should already be sorted by bit score, but just in case)
sort -k14,14nr ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats > ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats.sorted

mv ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats.sorted ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats

#echo ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats

#remove the blastdb files, fasta files for alignment
rm -f ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa.nhr
rm -f ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa.nin
rm -f ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa.nsq

rm -f ${out_path_mismatch_analysis}/${header_name}.chimp_coord.${seq_id}.fa
rm -f ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.fa

##file is blank, then no alignment appeared
if [[ ! -s ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats ]];
    then
        echo -e "empty_file\t${seq_id}"
else
    #num_lines=$(wc -l ${out_path}/${header_name}.human_coord.${seq_id}.blast.stats  | awk '{print $1}')
    #if [[ ${num_lines} -gt 1 ]];
        #then
            #echo -e "split_alignment\t${seq_id}"
    #fi
    ${python_dir}/python ${script_dir}/checkOligoFlankingMatchBlast.py ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats ${out_path_mismatch_analysis}/${header_name}.human_coord.${seq_id}.blast.stats.parsed

fi


