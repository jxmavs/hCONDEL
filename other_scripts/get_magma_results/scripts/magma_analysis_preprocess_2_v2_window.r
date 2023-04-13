#write out for GWAS overlap analyses
#see magma_analysis_preprocess_1_window.sh for generating hcondel_window_metadata

args = commandArgs(trailingOnly=TRUE)

all_genes_vec_file_path=args[1]
hcondel_window_metadata_file_path=args[2]
#size of window used to overlap genes
window_size=args[3]
additional_info=args[4]
output_file_path_with_header=args[5]
write_observed_overlap_genes=args[6]
seed_num=args[7]
write_random_overlap_genes=args[8]

####### read data #######

all_genes_vec=as.character(read.table(all_genes_vec_file_path, sep="\t")[,1])

hcondel_window_metadata=read.table(hcondel_window_metadata_file_path, header=F)

colnames(hcondel_window_metadata)=c("hCONDEL_ID", paste("hg38_protein_coding_window", window_size, "gene_name", sep="_"), paste("hg38_protein_coding_window", window_size, "gene_ensembl_id", sep="_"))

######### window output comparing against all hg38 protein coding genes ########

window_size=format(window_size, scientific = FALSE)
window_col=paste("hg38_protein_coding_window", window_size, "gene_ensembl_id", sep="_")

hcondel_window_metadata_filtered=hcondel_window_metadata

keep_ind=which(!is.na(hcondel_window_metadata_filtered[,window_col]))
overlap_genes=as.character(hcondel_window_metadata_filtered[keep_ind, window_col])
overlap_genes_list=strsplit(overlap_genes,"\\|")
overlap_genes_vec=unique(unlist(overlap_genes_list))
#print(length(overlap_genes_vec))

#remove genes which overlap in all_genes_vec
overlap_ind=which(all_genes_vec %in% overlap_genes_vec)
non_overlap_genes_vec=all_genes_vec[-overlap_ind]
#print(length(non_overlap_genes_vec))

#added 11/7/22, keep only overlap_genes that is in the list - ignore genes in which TSS does not map from hg38 to hg19
overlap_ind=which(overlap_genes_vec %in% all_genes_vec)
overlap_genes_vec=overlap_genes_vec[overlap_ind]

#separate out 
overlap_gene_out_df=data.frame(matrix(nrow=length(overlap_genes_vec)+length(non_overlap_genes_vec), ncol=2))
colnames(overlap_gene_out_df)=c("Gene_Id", "Group_Id")
overlap_gene_out_df[,1]=c(overlap_genes_vec, non_overlap_genes_vec)

overlap_gene_out_df[1:length(overlap_genes_vec),2]="Overlap"
overlap_gene_out_df[(length(overlap_genes_vec)+1):nrow(overlap_gene_out_df),2]="Not_Overlap"

if(write_observed_overlap_genes=="Y"){
	textOutFileName=paste(c(output_file_path_with_header, "window", window_size, additional_info, "hg38_protein_coding_gene_overlap_info.txt"), collapse="_")
	write.table(overlap_gene_out_df, textOutFileName, sep="\t", quote=F, row.names=F, col.names=F)
}

#write out random set matched based off overlap number
random_gene_out_df=overlap_gene_out_df
num_overlap=length(which(overlap_gene_out_df[,2]=="Overlap"))
random_gene_out_df[,2]="Not_Overlap"

if(write_random_overlap_genes=="Y"){
	#added 11/7/22, set.seed()
	set.seed(seed_num)	
	rand_ind=sample(nrow(random_gene_out_df), num_overlap)
	random_gene_out_df[rand_ind,2]="Overlap"

	textOutFileName=paste(c(output_file_path_with_header, "window", window_size, additional_info, "random_matched_hg38_protein_coding_gene_overlap_info.txt"), collapse="_")
	write.table(random_gene_out_df, textOutFileName, sep="\t", quote=F, row.names=F, col.names=F)
}
