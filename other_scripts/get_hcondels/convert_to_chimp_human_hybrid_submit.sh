
#this script creates the chimp human reference genomes for each coordinate
human_dir=/cluster_path/ape_project/deletions_project/reference_genome_files/human
filter_dir=/cluster_path/ape_project/deletions_project/filtering_new
chimp_human_hybrid_dir=/cluster_path/ape_project/deletions_project/chimp_human_hybrid

temp_dir=/cluster_path_temp/sgdp/temp
chimp_ref_file_names=/cluster_path/ape_project/deletions_project/chimp_human_hybrid/chimp_human_hybrid_genomes_files.txt
mkdir -p ${temp_dir}


rm -f ${temp_dir}/all_runs_convert_to_hybrid.txt
#submit jobs
while read -a array
do
    genome_file_id=${array[0]}
    #meta_file_name contains the chimp human hybrid coordinates in first part, then the chimp region in the second part    
    meta_file_name=${chimp_human_hybrid_dir}/${genome_file_id}_hybrid_out_coord.txt
    cons_overlap_del_file_name=${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_preconverted_${genome_file_id}.txt
    del_overlap_cons_file_name=${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_preconverted_${genome_file_id}.txt
    cons_partial_overlap_del_file_name=${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_preconverted_${genome_file_id}.txt    
    #first 3 coordinates are from original file?
    #last 6 columns contain new coordinates for where gap lies exactly, first 3 of that contain the conserved region, last 3 contains the exact gap region
    cons_overlap_del_converted_file_name=${chimp_human_hybrid_dir}/cleaned_conserved_overlapping_deletion_converted_${genome_file_id}.txt
    del_overlap_cons_converted_file_name=${chimp_human_hybrid_dir}/cleaned_deletion_overlapping_conserved_converted_${genome_file_id}.txt
    cons_partial_overlap_del_converted_file_name=${chimp_human_hybrid_dir}/cleaned_partial_conserved_overlapping_deletion_converted_${genome_file_id}.txt 

    echo "bash convert_to_chimp_human_hybrid_iter.sh ${meta_file_name}  ${cons_overlap_del_file_name} ${cons_overlap_del_converted_file_name} " >> ${temp_dir}/all_runs_convert_to_hybrid.txt
    echo "bash convert_to_chimp_human_hybrid_iter.sh ${meta_file_name}  ${del_overlap_cons_file_name} ${del_overlap_cons_converted_file_name} " >> ${temp_dir}/all_runs_convert_to_hybrid.txt
    echo "bash convert_to_chimp_human_hybrid_iter.sh ${meta_file_name}  ${cons_partial_overlap_del_file_name} ${cons_partial_overlap_del_converted_file_name} " >> ${temp_dir}/all_runs_convert_to_hybrid.txt
 
done<${chimp_ref_file_names}

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project 
totalMem=3G
totalNumTasks=$(wc -l ${temp_dir}/all_runs_convert_to_hybrid.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEval.sh > ${temp_dir}/submit_convert_to_hybrid.sh
qsub ${temp_dir}/submit_convert_to_hybrid.sh ${temp_dir}/all_runs_convert_to_hybrid.txt

