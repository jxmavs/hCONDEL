
ID=$1
seqWithoutEdit=$2
seqWithEdit=$3
guideRNA=$4
cleavageOffset=$5
inputFilePath=$6
outputFilePath=$7

#CRISPResso -r1 ${inputFilePath}/${ID}_R1.fastq.gz -r2 ${inputFilePath}/${ID}_R2.fastq.gz -a ${seqWithoutEdit} -e ${seqWithEdit} -g ${guideRNA} --cleavage_offset ${cleavageOffset} --hdr_perfect_alignment_threshold 98 --output_folder ${outputFilePath}
#python -c "import sys; print(sys.executable)"
#CRISPResso2

#CRISPResso --fastq_r1 ${inputFilePath}/${ID}_R1.fastq.gz --fastq_r2 ${inputFilePath}/${ID}_R2.fastq.gz --amplicon_seq ${seqWithoutEdit} --guide_seq ${guideRNA} --expected_hdr_amplicon_seq ${seqWithEdit} --output_folder ${outputFilePath}

#CRISPResso --fastq_r1 ${inputFilePath}/${ID}_R1.fastq.gz --fastq_r2 ${inputFilePath}/${ID}_R2.fastq.gz --amplicon_seq ${seqWithoutEdit} --guide_seq ${guideRNA} --expected_hdr_amplicon_seq ${seqWithEdit} --quantification_window_size 0 --quantification_window_center ${cleavageOffset}  --min_average_read_quality 30 --min_single_bp_quality 10 --min_paired_end_reads_overlap 7 --output_folder ${outputFilePath}

#CRISPResso --fastq_r1 ${inputFilePath}/${ID}_R1.fastq.gz --fastq_r2 ${inputFilePath}/${ID}_R2.fastq.gz --amplicon_seq ${seqWithoutEdit} --guide_seq ${guideRNA} --expected_hdr_amplicon_seq ${seqWithEdit} --quantification_window_size 20 --quantification_window_center ${cleavageOffset}  --min_average_read_quality 20 --min_single_bp_quality 10 --min_paired_end_reads_overlap 7 --output_folder ${outputFilePath}

CRISPResso --fastq_r1 ${inputFilePath}/${ID}_R1.fastq.gz --fastq_r2 ${inputFilePath}/${ID}_R2.fastq.gz --amplicon_seq ${seqWithoutEdit} --guide_seq ${guideRNA} --expected_hdr_amplicon_seq ${seqWithEdit} --quantification_window_size 20 --quantification_window_center ${cleavageOffset}  --min_average_read_quality 30 --min_single_bp_quality 10 --min_paired_end_reads_overlap 7 --output_folder ${outputFilePath}
