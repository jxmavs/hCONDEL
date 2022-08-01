## Main Prereqs
* python
* bwa
* R
* UCSC tools

* Main R packages:
  * Seurat
  * ggplot2
  * DESeq2

## Analyze processed data
See analyze_processed_data/hCONDEL_processed_data_analysis_script.r for generating figures relating to the DESeq2 processed data.
See analyze_processed_data/single_cell_analyses_del_process.r for generating figures relating to the LOXL2-edited cells.

## Generating multiple sequence alignment
See other_scripts/generate_multiz_alignment for scripts used to generate the 11-way multiple sequence alignment in the paper. multiple_sequence_master_script.sh is the main script. The file all.mod in the directory is the neutral parameter background file for running phastCons.

## Getting hCONDELs
See other_scripts/get_hcondels for scripts used to obtain hCONDELs from the 11-way multiple sequence alignment/pairwise alignment between human (hg38) and chimpanzee (panTro4). finalFilteringProtocol_v4.sh, followed by deletion_protocol_continued.sh are the main scripts.

## Preprocess read data
See other_scripts/get_counts for scripts used to process read data to generate count files. rep_barcode_associate_full_dataset.sh is the main script.

## Generate processed data
After generating the initial count data, to generate the processed data (DESeq2 fold changes/p-values), see other_scripts/get_deseq2_data. 

## Contact
Feel free to contact jxue@broadinstitute.org for any other inquires/questions.
