#uncomment each of the sections to get results for each analysis

temp_dir=/cluster_path/ape_project/deletions_project/hdr_editting_exp/temp

############## (LOXL2 edit proportions after HDR in SK-N-SHs) 2_19_20 #####################
#input_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/files/2_19_20_run
#output_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/output/2_19_20_run/orig_quant_window
#mkdir -p ${output_file_path}
#CRISPResso_ref_file=${input_file_path}/hdr_editting_exp_CRISPResso_input_2_19_20.txt
#run_name=CRISPResso_2_19_20_del
#################################################

############## (HCR-FlowFish on LOXL2 HDR-edited SK-N-SHs, rep 1) 8_4_20 #####################
#input_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/files/8_4_20_run
#output_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/output/8_4_20_run/orig_quant_window
#mkdir -p ${output_file_path}
#CRISPResso_ref_file=${input_file_path}/hdr_editting_exp_CRISPResso_input_8_4_20.txt
#run_name=CRISPResso_8_4_20_del
#################################################

############## (HCR-FlowFish on LOXL2 HDR-edited SK-N-SHs, rep 2) 9_8_20 #####################
#input_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/files/9_8_20_run
#output_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/output/9_8_20_run/orig_quant_window
#mkdir -p ${output_file_path}
#CRISPResso_ref_file=${input_file_path}/hdr_editting_exp_CRISPResso_input_9_8_20.txt
#run_name=CRISPResso_9_8_20_del
#################################################

############## (PPP2CA promoter CRISPR-mutagenesis in SK-N-SHs, rep 1-3) 11_24_22 #####################
#input_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/files/11_24_22_run
#output_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/output/11_24_22_run/orig_quant_window
#output_file_path=/cluster_path/ape_project/deletions_project/hdr_editting_exp/output/11_24_22_run/orig_quant_window_q20
#mkdir -p ${output_file_path}
#CRISPResso_ref_file=${input_file_path}/hdr_editting_exp_CRISPResso_input_11_24_22.txt
#run_name=CRISPResso_11_24_22_del
#################################################

script_dir=/cluster_path/ape_project/deletions_project/hdr_editting_exp/scripts
#remove header
awk 'NR>1{print}' ${CRISPResso_ref_file} > ${temp_dir}/temp_ref_file.txt


rm -f ${temp_dir}/${run_name}_all_runs.txt
#submit jobs
while read -a array
do
    ID=${array[0]}
    seqWithoutEdit=${array[1]}
    seqWithEdit=${array[2]}
    guideRNA=${array[3]}
    cleavageOffset=${array[4]}
    echo "bash ${script_dir}/run_CRISPResso_iter.sh ${ID} ${seqWithoutEdit} ${seqWithEdit} ${guideRNA} ${cleavageOffset} ${input_file_path} ${output_file_path}" >> ${temp_dir}/${run_name}_all_runs.txt
   
done<${temp_dir}/temp_ref_file.txt

#submit all jobs
del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${temp_dir}/${run_name}_all_runs.txt  | awk '{print $1}')
totalTime=10:00:00
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l h_vmem=totalMem/#$ -l h_vmem=${totalMem}/g; s/#$ -l h_rt=totalTime/#$ -l h_rt=${totalTime}/g" ${del_dir}/runTasksEvalLDRevised.sh > ${temp_dir}/${run_name}_submit_all.sh
rm -f ${temp_dir}/${run_name}_submit_all.log
qsub -o ${temp_dir}/${run_name}_submit_all.log ${temp_dir}/${run_name}_submit_all.sh ${temp_dir}/${run_name}_all_runs.txt
