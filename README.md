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
See other_scripts/generate_multiz_alignment for scripts used to generate the 11-way multiple sequence alignment in the paper. multiple_sequence_master_script.sh is the main script.

## Getting hCONDELs
See other_scripts/get_hcondels for scripts used to obtain hCONDELs from the 11-way multiple sequence alignment/pairwise alignment between human (hg38) and chimpanzee (panTro4). finalFilteringProtocol_v4.sh, followed by deletion_protocol_continued.sh are the main scripts.

See denylist_commands_with_hg18_check_added_panTro6.sh for further filtering hCONDELs to the 10,032 hCONDELs used in the paper. This filtering script includes overlapping the hCONDELs with the GAGP dataset and checking that hCONDELs can be mapped to panTro5 (see materials and methods in the paper).

## Generate the barcode table
See other_scripts/build_barcode_table for scripts used to generate the barcode table, "enhancer_barcode_1_22_19_merged_cutAMTBegin_16_adapterSeqTrimEnd_trimmed_bowtie2_barcode_assoc_table.txt" which is used in "other_scripts/get_counts/rep_barcode_associate_full_dataset.sh"

other_scripts/build_barcode_table/plasmid_seq_analyze_cleaned.sh is the main script.

## Preprocess read data
See other_scripts/get_counts for scripts used to process read data to generate count files. rep_barcode_associate_full_dataset.sh is the main script.

## Generate processed data
After generating the initial count data, to generate the processed data (DESeq2 fold changes/p-values), see other_scripts/get_deseq2_data/run_deseq2.r 

## hCONDEL TF motif alteration scores
See other_scripts/get_tf_alteration_scores_and_enrichment for scripts for calculating the TF alteration scores.

other_scripts/get_tf_alteration_scores_and_enrichment/scripts/get_tf_alteration_scores_JASPAR_2020.sh is the main script.

## TF motif clustering
See other_scripts/get_tf_alteration_scores_and_enrichment/tf_clustering for scripts for generating the TF clusters for the JASPAR 2020 motifs. It is adapted from Vierstra et al. 2020.

other_scripts/get_tf_alteration_scores_and_enrichment/tf_clustering/scripts/Workflow_v2.1beta-human.modified.ipynb is the main script.

## hCONDEL TF alteration enrichments
See other_scripts/get_tf_alteration_scores_and_enrichment for scripts for generating the hCONDEL TF motif perturbation enrichments.

other_scripts/get_tf_alteration_scores_and_enrichment/scripts/tf_enrichment_analysis_submission.sh is the main script.

## hCONDEL GWAS enrichments
See other_scripts/get_magma_results for analyses scripts for calculating the GWAS enrichments.

other_scripts/get_magma_results/scripts/magma_analysis_v3_window.sh is the main script that can generate the GWAS enrichments for the displayed neurological traits in Fig. 1G.

other_scripts/get_magma_results/scripts/magma_analysis_ukbb_v2_window.sh is the main script that can generate the hCONDEL GWAS enrichments from 4,178 traits from the UKBB

## hCONDEL genomic annotation and GTEx enrichments
See other_scripts/get_genomic_annot_GTEx_enrichments for scripts for calculating enrichments for hCONDEL age, overlap with different repeat classes, coding, cis-regulatory element, and GTEx RNA-seq.

other_scripts/get_genomic_annot_GTEx_enrichments/scripts/hcondel_other_enrichment.sh is the main script.

## hCONDEL H3K27ac enrichments
See other_scripts/get_H3K27ac_roadmap_enrichments for scripts to calculate hCONDEL enrichments on H3K27ac peaks from the Roadmap Epigenomics Project.

other_scripts/get_H3K27ac_roadmap_enrichments/scripts/analyze_roadmap_data_v5_submission.sh is the main script.

## CRISPResso analyses
See other_scripts/crispresso_analyses for scripts for getting the allelic proportions from the LOXL2 HDR edit, the LOXL2 HCR-FlowFish, and the PPP2CA promoter mutagenesis experiments.

other_scripts/crispresso_analyses/run_CRISPResso_submit_submission.sh is the main script.

## LOXL2 single-cell preprocess
See other_scripts/loxl2_single_cell_preprocess for scripts for processing the single cell RNA-Seq data from LOXL2-edited SK-N-SHs. Scripts for processing the enriched LOXL2-edited locus reads for classifying SK-N-SHs with and without the chimpanzee-inserted base is also found here.

other_scripts/loxl2_single_cell_preprocess/cell_ranger/cell_ranger_commands_submission.sh is the main script for processing the single cell RNA-Seq data from LOXL2-edited SK-N-SHs following the 10X pipeline.

other_scripts/loxl2_single_cell_preprocess/got/got_commands_submission.sh is the main script for processing the enriched LOXL2-edited locus reads using the GoT pipeline.

## Contact
Feel free to contact jxue@broadinstitute.org for any other inquires/questions.
