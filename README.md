## Prereqs
* python
* bwa
* R

* R packages:
  * ggplot2
  * DESeq2

## Analyze processed data
See scripts/analyze_processed_data/hCONDEL_processed_data_analysis_script.r for generating figures relating to the DESeq2 processed data.
See scripts/analyze_processed_data/single_cell_analyses_del_process.r for generating figures relating to the LOXL2-edited cells.

## Generating Multiz alignment/getting hCONDELs
See scripts/generate_multiz_alignment for scripts used to generate the 11-way multiz alignment in the paper

## Getting hCONDELs
See scripts/get_hcondels for scripts used to obtain hCONDELs from the 11-way multiple sequence alignment/pairwise alignment between human (hg38) and chimpanzee (panTro4)

## Preprocess read data
See scripts/get_counts for scripts used to process read data to generate count files

## Generate processed data
After generating the initial count data, to generate the processed data (DESeq2 fold changes/p-values), see scripts/files in scripts/get_deseq2_data 

