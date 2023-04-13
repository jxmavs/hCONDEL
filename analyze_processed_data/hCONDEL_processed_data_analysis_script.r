library(ggplot2)
suppressMessages(library(DESeq2))
library(ggsignif)
library(ggrepel)
library(pheatmap)

#******************** read in original dataset (preprocessing, run before any analyses) ********************#

input_file_path="/home_path/NovaSeq_Analysis/Analysis_Output"
output_file_path="/home_path/NovaSeq_Analysis/Analysis_Output"
out_prefix="hCONDEL"
output_file_path_with_header=paste(c(output_file_path, out_prefix), collapse="/")

#read data, from Supplementary Table 1

#read in metainfo

hcondel_mpra_metadata_file_name="novaseq_12_15_18_hCONDEL_cleaned_metatable_submission.txt"
hcondel_mpra_metadata_file_path=paste(input_file_path, hcondel_mpra_metadata_file_name, sep="/")
hcondel_mpra_metadata=read.table(hcondel_mpra_metadata_file_path, header=T)

#read in deseq2 results

hcondel_mpra_deseq2_results_file_name="novaseq_12_15_18_deseq_results_submission.txt"
hcondel_mpra_deseq2_results_file_path=paste(input_file_path, hcondel_mpra_deseq2_results_file_name, sep="/")
hcondel_mpra_deseq2_results=read.table(hcondel_mpra_deseq2_results_file_path, header=T)

#read in the plasmid count file 

hcondel_mpra_plasmid_count_table_file_name="novaseq_12_15_18_count_data_submission.txt"
hcondel_mpra_plasmid_count_table_file_path=paste(input_file_path, hcondel_mpra_plasmid_count_table_file_name, sep="/")
hcondel_mpra_plasmid_count_table=read.table(hcondel_mpra_plasmid_count_table_file_path, header=T)

###### get average of plasmid counts, use it for filtering later #######

human_plasmid_count_cols=paste("Plasmid", c("R1", "R2", "R3", "R4", "R5"), "Human", sep="_")
chimp_plasmid_count_cols=paste("Plasmid", c("R1", "R2", "R3", "R4", "R5"), "Chimp", sep="_")

plasmid_avg_human=apply(hcondel_mpra_plasmid_count_table[, human_plasmid_count_cols], 1, function(x){mean(x, na.rm=T)})
plasmid_avg_chimp=apply(hcondel_mpra_plasmid_count_table[, chimp_plasmid_count_cols], 1, function(x){mean(x, na.rm=T)})

#make aggregate dataset for analysis
hcondel_full_dataset=cbind(hcondel_mpra_metadata, hcondel_mpra_deseq2_results[,-1], plasmid_avg_human, plasmid_avg_chimp)

#find low plasmid count indeces for filtering later

plasmid_cutoff_map=list()
plasmid_cutoff_map[["HEK293"]]=20
plasmid_cutoff_map[["HEPG2"]]=20
plasmid_cutoff_map[["K562"]]=20
plasmid_cutoff_map[["GM12878"]]=20
plasmid_cutoff_map[["SKNSH"]]=20
plasmid_cutoff_map[["NPC"]]=60

low_plasmid_count_indeces_human_all=list()
low_plasmid_count_indeces_chimp_all=list()
low_plasmid_count_indeces_human_and_chimp_all=list()

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")
for(ind in 1:length(cell_names)){

	cell_name=cell_names[ind]
	plasmid_cutoff_val=plasmid_cutoff_map[[cell_name]]
	
	low_plasmid_count_indeces_human=which( hcondel_full_dataset[, "plasmid_avg_human"]<plasmid_cutoff_val )
	low_plasmid_count_indeces_human_all[[cell_name]]=low_plasmid_count_indeces_human

	low_plasmid_count_indeces_chimp=which( hcondel_full_dataset[, "plasmid_avg_chimp"]<plasmid_cutoff_val )
	low_plasmid_count_indeces_chimp_all[[cell_name]]=low_plasmid_count_indeces_chimp	
	
	low_plasmid_count_indeces_human_and_chimp=which( (hcondel_full_dataset[, "plasmid_avg_human"]<plasmid_cutoff_val) | (hcondel_full_dataset[, "plasmid_avg_chimp"]<plasmid_cutoff_val) )
	low_plasmid_count_indeces_human_and_chimp_all[[cell_name]]=low_plasmid_count_indeces_human_and_chimp
	
}

#******************** GREAT bar plot (Figure 1F) ********************#

#from the original GREAT output, "#" is removed in the header line and saved, this allows the header to be read
great_file_name="/home_path/hCONDEL_del_greatExportAll_7_11_22.tsv"
great_stats=read.table(great_file_name, sep="\t", comment.char="#", header=T)

great_stats=great_stats[1:15,]

great_stats[,"logFDRQValue"]=-log(great_stats[, "BinomFdrQ"], base=10)

#sort by logFDRQValue in plotting
sortedInd=order(great_stats[,"logFDRQValue"])
great_stats=great_stats[sortedInd,]
great_stats[,"Desc"]=factor(great_stats[,"Desc"], as.character(great_stats[,"Desc"]))

p=ggplot(great_stats, aes_string(x = "Desc", y = "logFDRQValue")) +geom_bar(stat = "identity", fill = "yellowgreen") +coord_flip()+ylab("Log Q Value")+xlab("")+theme(axis.text.x = element_text(size=60), axis.title.x = element_text(size=50), axis.text.y = element_text(size=40), axis.title.y = element_text(size=40), panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotName="great_deletion_comp_whole_genome_bar_chart"
plotOutFileName=paste(c(output_file_path_with_header,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
ggsave(filename=plotOutFileName, plot=p, width=35, height=20)

#*********************************** fimo tf difference analysis (Fig. 2C and fig. S6B,C) ********************************************************************#

##### parameters to set #####

#analysis type, set to either "MPRA_FIMO_correlation" (used for fig. S6B), "MPRA_FIMO_correlation_pos_active_in_at_least_one_context" (used to create pie chart for Fig. 2C), or "MPRA_FIMO_correlation_no_cons" (used to get the correlation sats for fig. S6C)
analysis_type="MPRA_FIMO_correlation"
#analysis_type="MPRA_FIMO_correlation_pos_active_in_at_least_one_context"
#analysis_type="MPRA_FIMO_correlation_no_cons"
#include TF labels or no, set to false if don't want them

include_tf_labels=TRUE

#############################

#first entry is activity padj threshold, second is skew padj threshold
#set padj threshold to 0.2 for all cell types except NPCs, which is set to 0.05
padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.2)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.2)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.2)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.2)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.2)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.1)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.1)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.1)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.1)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.1)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

save_indiv_plots=FALSE

#filters to remove hCONDELs with an excess of mismatches between the human and chimp sequence
num_bp_unaligned_cutoff=quantile(hcondel_full_dataset[,"num_bp_unaligned"], probs = seq(0, 1, 0.1))[10]
num_bp_mismatches_cutoff=quantile(hcondel_full_dataset[,"num_bp_mismatches"], probs = seq(0, 1, 0.1))[9]

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

tf_header_name="predicted_JASPARv2020_tf"
tf_analysis_subset_table_combined=data.frame()

##############################

#function to obtain format of correlation and agree/disagree
get_cor_stats = function(df, y_name, x_name){
	
	num_agree=length(which(sign(df[,y_name])==sign(df[,x_name])))
	num_disagree=length(which(sign(df[,y_name])!=sign(df[,x_name])))
	
	cor_test_results=cor.test(df[,y_name], df[,x_name], method="pearson")
	r=format(cor_test_results$estimate, digits=2)
	p=format(cor_test_results$p.value, digits=2)
	
	#binomial test assuming 50% is true
	agree_pval=format(binom.test(num_agree, num_agree+num_disagree, p=0.5)$p.value, digits=2)
	
	out_str=sprintf("r=%s\np-val=%s\nagree:%s\ndisagree:%s\nagreement p-val:%s", r, p, num_agree, num_disagree, agree_pval)
	
	out_list=list()
	out_list[[1]]=out_str
	out_list[[2]]=num_agree
	out_list[[3]]=num_disagree
	
	return(out_list)
}


for( i in 1:length(cell_names)){

	cell_name=cell_names[i]
	
	low_plasmid_count_indeces_human_and_chimp=low_plasmid_count_indeces_human_and_chimp_all[[cell_name]]
	
	#remove low count indeces and hCONDELs with an excess of mismatches between the human and chimp sequence
	hcondel_full_dataset_filtered=hcondel_full_dataset[-low_plasmid_count_indeces_human_and_chimp, ]
	hcondel_full_dataset_filtered=hcondel_full_dataset_filtered[which( (hcondel_full_dataset_filtered[,"num_bp_unaligned"]<num_bp_unaligned_cutoff) & (hcondel_full_dataset_filtered[,"num_bp_mismatches"]<num_bp_mismatches_cutoff)),]

	padj_activity_cutoff=padj_val_sig_cutoff_list[[cell_name]][1]
	padj_skew_cutoff=padj_val_sig_cutoff_list[[cell_name]][2]
	
	lfc_chimp_col=paste("log2FoldChange_Chimp", cell_name, sep="_")
	lfc_human_col=paste("log2FoldChange_Human", cell_name, sep="_")
	
	lfc_skew_col=paste("log2FoldChange_Skew", cell_name, sep="_")
	
	padj_chimp_col=paste("padj_Chimp", cell_name, sep="_")
	padj_human_col=paste("padj_Human", cell_name, sep="_")
	
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")
	
	tf_header_name_cell=paste(cell_name, tf_header_name, sep="_")
	
	tf_name=paste(tf_header_name_cell, "name", sep="_")
	
	hg38_panTro4_tf_diff_col=paste(tf_header_name_cell, "hg38_panTro4_diff", sep="_")
	hg38_rheMac8_tf_diff_col=paste(tf_header_name_cell, "hg38_rheMac8_diff", sep="_")

	hg38_tf_other_stats_col=paste(paste(tf_header_name_cell, "hg38", sep="_"), c("score", "pval"), sep="_")
	panTro4_tf_other_stats_col=paste(paste(tf_header_name_cell, "panTro4", sep="_"), c("score", "pval"), sep="_")
	rheMac8_tf_other_stats_col=paste(paste(tf_header_name_cell, "rheMac8", sep="_"), c("score", "pval"), sep="_")

	############################################
	
	if(analysis_type=="MPRA_FIMO_correlation"){
		#filters based off of: zoonomia phyloP scores, strong repressive activity in human and chimp context, phastCons conserved block score, require sign of change of TF difference to be the same between between panTro4 and rheMac8
		keep_ind=which( !(is.na(hcondel_full_dataset_filtered[,tf_name])) & (hcondel_full_dataset_filtered[, lfc_chimp_col]>-0.5) & (hcondel_full_dataset_filtered[, lfc_human_col]>-0.5) & (hcondel_full_dataset_filtered[, padj_skew_col]<padj_skew_cutoff) & (hcondel_full_dataset_filtered[,"del_phyloP_max_score"]>1) & (hcondel_full_dataset_filtered[,"cons_phastCons_score"]>50) & (sign(hcondel_full_dataset_filtered[, hg38_panTro4_tf_diff_col])==sign(hcondel_full_dataset_filtered[, hg38_rheMac8_tf_diff_col ])) )
	} else if(analysis_type=="MPRA_FIMO_correlation_pos_active_in_at_least_one_context") {
		#For tf alteration category estimates in Fig. 2C, zoonomia phyloP scores, positive enhancer activity in either the human or chimp context, phastCons conserved block score, require sign of change of TF difference to be the same between between panTro4 and rheMac8
		keep_ind=which( !(is.na(hcondel_full_dataset_filtered[,tf_name])) & ( ((hcondel_full_dataset_filtered[, lfc_chimp_col]>0) & (hcondel_full_dataset_filtered[, padj_chimp_col]<padj_activity_cutoff)) | ((hcondel_full_dataset_filtered[, lfc_human_col]>0) & (hcondel_full_dataset_filtered[, padj_human_col]<padj_activity_cutoff)) ) & (hcondel_full_dataset_filtered[, padj_skew_col]<padj_skew_cutoff) & (hcondel_full_dataset_filtered[,"del_phyloP_max_score"]>1) & (hcondel_full_dataset_filtered[,"cons_phastCons_score"]>50) & (sign(hcondel_full_dataset_filtered[, hg38_panTro4_tf_diff_col])==sign(hcondel_full_dataset_filtered[, hg38_rheMac8_tf_diff_col ])) )
	} else if(analysis_type=="MPRA_FIMO_correlation_no_cons"){
		#same as MPRA_FIMO_correlation, except with the zoonomia phyloP score filter removed
		keep_ind=which( !(is.na(hcondel_full_dataset_filtered[,tf_name])) & (hcondel_full_dataset_filtered[, lfc_chimp_col]>-0.5) & (hcondel_full_dataset_filtered[, lfc_human_col]>-0.5) & (hcondel_full_dataset_filtered[, padj_skew_col]<padj_skew_cutoff) & (hcondel_full_dataset_filtered[,"cons_phastCons_score"]>50) & (sign(hcondel_full_dataset_filtered[, hg38_panTro4_tf_diff_col])==sign(hcondel_full_dataset_filtered[, hg38_rheMac8_tf_diff_col ])) )
	} else{
		print("Invalid Analysis Type!")
		break
	}
	
	###############################################################

	tf_analysis_subset_table=hcondel_full_dataset_filtered[keep_ind,]
	tf_analysis_subset_table[, tf_name]=factor(as.character(tf_analysis_subset_table[, tf_name]))

	### add to combined cell type table, keep relevant columns ###
	keep_cols=c("hCONDEL_ID", lfc_skew_col, padj_skew_col, tf_name, hg38_panTro4_tf_diff_col, hg38_rheMac8_tf_diff_col, lfc_chimp_col, padj_chimp_col, lfc_human_col, padj_human_col, hg38_tf_other_stats_col, panTro4_tf_other_stats_col, rheMac8_tf_other_stats_col)
	tf_analysis_subset_table_iter=tf_analysis_subset_table[,keep_cols]
	
	#remove cell name in columns of combined table
	colnames(tf_analysis_subset_table_iter)=gsub(paste(cell_name, "_", sep=""), "", colnames(tf_analysis_subset_table_iter))
	colnames(tf_analysis_subset_table_iter)=gsub(paste("_", cell_name, sep=""), "", colnames(tf_analysis_subset_table_iter))
	tf_analysis_subset_table_combined=rbind(tf_analysis_subset_table_combined, tf_analysis_subset_table_iter)

	#if want to save individual plot per cell type
	if(save_indiv_plots & (nrow(tf_analysis_subset_table)>1)){
		x_col_name=hg38_panTro4_tf_diff_col
		y_col_name=paste("log2FoldChange_Skew", cell_name, sep="_")
	
		plot_name=paste(cell_name, analysis_type, "human-chimp_skew_vs_FIMO_binding_score_diff", sep="_")
		plot_output_file_path=paste(output_file_path_with_header, paste(plot_name, ".pdf", sep=""), sep="_")
	
		get_cor_stats_output=get_cor_stats(tf_analysis_subset_table, y_col_name, x_col_name)
		cor_stats_str=get_cor_stats_output[[1]]
		print(cor_stats_str)
	
		p_hg38_panTro4=ggplot(tf_analysis_subset_table, aes_string(x=x_col_name, y=y_col_name, label=tf_name))+geom_point(alpha=0.7, size=4)+
		xlab("Human/Chimp Max tf Binding Score Difference")+geom_text_repel(size=4, max.overlaps = Inf)+ylab("log(Human/Chimp) Skew")+
		geom_vline(xintercept = 0, linetype="dashed")+geom_hline(yintercept = 0, linetype="dashed")+
		stat_smooth(method="lm", se=FALSE)+annotate("text", x = 15, y = -1.5, label=cor_stats_str, size=4)+
		theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), 
		axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))

		ggsave(filename=plot_output_file_path, plot=p_hg38_panTro4 )
	}
}

#################### make combined plot #########################
#run immediately after making the individual plots before

#if there are duplicate hCONDEL IDs, keep the one with the most significant, as measured by padj skew

tf_name_combined=paste(tf_header_name, "name", sep="_")

hg38_panTro4_tf_diff_col_combined=paste(tf_header_name, "hg38_panTro4_diff", sep="_")
hg38_rheMac8_tf_diff_col_combined=paste(tf_header_name, "hg38_rheMac8_diff", sep="_")

hg38_tf_other_stats_col_combined=paste(paste(tf_header_name, "hg38", sep="_"), c("score", "pval"), sep="_")
panTro4_tf_other_stats_col_combined=paste(paste(tf_header_name, "panTro4", sep="_"), c("score", "pval"), sep="_")
rheMac8_tf_other_stats_col_combined=paste(paste(tf_header_name, "rheMac8", sep="_"), c("score", "pval"), sep="_")

unique_ids=names(table(as.character(tf_analysis_subset_table_combined[,"hCONDEL_ID"])))

tf_analysis_subset_table_combined_unique=data.frame(matrix(nrow=length(unique_ids), ncol=ncol(tf_analysis_subset_table_combined)))
colnames(tf_analysis_subset_table_combined_unique)=colnames(tf_analysis_subset_table_combined)

for(ind in 1:length(unique_ids)){
	unique_id_iter=unique_ids[ind]
	subset_ind=which(tf_analysis_subset_table_combined[,"hCONDEL_ID"]==unique_id_iter)
	
	most_sig_ind=which.min(tf_analysis_subset_table_combined[subset_ind, "padj_Skew"])
	
	tf_analysis_subset_table_combined_unique[ind, c("hCONDEL_ID", tf_name_combined)]=as.character(unlist(tf_analysis_subset_table_combined[subset_ind,][most_sig_ind, c("hCONDEL_ID", tf_name_combined)]))
	
	numeric_cols=c("log2FoldChange_Skew", "padj_Skew", hg38_panTro4_tf_diff_col_combined, hg38_rheMac8_tf_diff_col_combined, hg38_tf_other_stats_col_combined, panTro4_tf_other_stats_col_combined, rheMac8_tf_other_stats_col_combined)
	tf_analysis_subset_table_combined_unique[ind ,numeric_cols]=tf_analysis_subset_table_combined[subset_ind,][most_sig_ind, numeric_cols]

}

### hg38 compared with panTro4 ###

x_col_name=hg38_panTro4_tf_diff_col_combined
y_col_name="log2FoldChange_Skew"

get_cor_stats_output=get_cor_stats(tf_analysis_subset_table_combined_unique, y_col_name, x_col_name)
cor_stats_str=get_cor_stats_output[[1]]
print(cor_stats_str)

if(include_tf_labels){
	#all plot
	label_info="with_labels"
	p_hg38_panTro4=ggplot(tf_analysis_subset_table_combined_unique, aes_string(x=x_col_name, y=y_col_name, label=tf_name_combined))+geom_point(alpha=0.7, size=4)+
		xlab("Human/Chimp Max tf Binding Score Difference")+geom_text_repel(size=4, max.overlaps = Inf)+ylab("log(Human/Chimp) Skew")+
		geom_vline(xintercept = 0, linetype="dashed")+geom_hline(yintercept = 0, linetype="dashed")+
		stat_smooth(method="lm", se=FALSE)+annotate("text", x = 15, y = -2.25, label=cor_stats_str, size=4)+
		theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), 
		axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
} else{
	#no labels
	label_info="no_labels"
	p_hg38_panTro4=ggplot(tf_analysis_subset_table_combined_unique, aes_string(x=x_col_name, y=y_col_name))+geom_point(alpha=0.7, size=4)+
		xlab("Human/Chimp Max tf Binding Score Difference")+ylab("log(Human/Chimp) Skew")+
		geom_vline(xintercept = 0, linetype="dashed")+geom_hline(yintercept = 0, linetype="dashed")+
		stat_smooth(method="lm", se=FALSE)+annotate("text", x = 15, y = -2.25, label=cor_stats_str, size=4)+
		theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), 
		axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
}

plot_name=paste(analysis_type, "human-chimp_skew_vs_FIMO_binding_score_diff_combined", label_info, sep="_")
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p_hg38_panTro4 )


####### make pie chart for activators/repressors ########

num_create_or_strengthen_activator=length(which((tf_analysis_subset_table_combined_unique[,"log2FoldChange_Skew"]>0) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_panTro4_diff"]>0)))
num_destroy_or_weaken_activator=length(which((tf_analysis_subset_table_combined_unique[,"log2FoldChange_Skew"]<0) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_panTro4_diff"]<0)))
num_create_or_strengthen_repressor=length(which((tf_analysis_subset_table_combined_unique[,"log2FoldChange_Skew"]<0) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_panTro4_diff"]>0)))
num_destroy_or_weaken_repressor=length(which((tf_analysis_subset_table_combined_unique[,"log2FoldChange_Skew"]>0) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_panTro4_diff"]<0)))
total_num=num_create_or_strengthen_activator+num_destroy_or_weaken_activator+num_create_or_strengthen_repressor+num_destroy_or_weaken_repressor

print(total_num==nrow(tf_analysis_subset_table_combined_unique))

prop_create_or_strengthen_activator=num_create_or_strengthen_activator/total_num
prop_destroy_or_weaken_activator=num_destroy_or_weaken_activator/total_num
prop_create_or_strengthen_repressor=num_create_or_strengthen_repressor/total_num
prop_destroy_or_weaken_repressor=num_destroy_or_weaken_repressor/total_num

categories_revised_names=c("Create_or_Strengthen_Activator", "Destroy_or_Weaken_Activator", "Create_or_Strengthen_Repressor", "Destroy_or_Weaken_Repressor")

pie_chart_df_prop=data.frame(matrix(NA, nrow=length(categories_revised_names), ncol=3))
colnames(pie_chart_df_prop)=c("class", "prop", "y.pos")
pie_chart_df_prop[,"prop"]=round(c(prop_create_or_strengthen_activator, prop_destroy_or_weaken_activator, prop_create_or_strengthen_repressor, prop_destroy_or_weaken_repressor), 3)
pie_chart_df_prop[,"class"]=categories_revised_names

#sort by prop
sorted_ind=order(pie_chart_df_prop[,"prop"], decreasing=TRUE)

pie_chart_df_prop=pie_chart_df_prop[sorted_ind,]

mycols <- c("green3", "gold2", "tomato1", "gray70")

pie_chart_df_prop$class <- factor(pie_chart_df_prop$class, pie_chart_df_prop$class)

pie_chart_df_prop[,"y.pos"]=1-(cumsum(c(0, pie_chart_df_prop[,"prop"])) + c(pie_chart_df_prop[,"prop"] / 2, .01))[1:nrow(pie_chart_df_prop)]

p=ggplot(pie_chart_df_prop, aes(x = "", y = prop, fill = class)) +geom_bar(width = 1, stat = "identity", color = "white")+coord_polar("y", start = 0)+
	geom_label_repel(aes(y=y.pos, label = prop), size=40, show.legend = F, nudge_x = 2, nudge_y = 20)+guides(fill = guide_legend(title = ""), box.padding = unit(5, "lines"))+scale_fill_manual(values = mycols)+
	theme(legend.text=element_text(size=55, color="black"), legend.title=element_text(size=25), legend.key.size = unit(3, "cm"),  
	axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
	panel.border = element_blank(), panel.grid=element_blank(), panel.background = element_blank())

plot_name=paste(analysis_type, "activator_repressor", "pie_chart", sep="_")
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p, width=25, height=25)

####### make pie chart for creating novel motifs ########

p_val_cutoff=0.0001

num_create_motif=length(which((tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_pval"]<p_val_cutoff) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_panTro4_pval"]>=p_val_cutoff)))
num_destroy_motif=length(which((tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_pval"]>=p_val_cutoff) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_panTro4_pval"]<p_val_cutoff)))
num_modulate_motif=length(which((tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_hg38_pval"]<p_val_cutoff) & (tf_analysis_subset_table_combined_unique[,"predicted_JASPARv2020_tf_panTro4_pval"]<p_val_cutoff)))

total_num=num_create_motif+num_destroy_motif+num_modulate_motif

prop_create_motif=num_create_motif/total_num
prop_destroy_motif=num_destroy_motif/total_num
prop_modulate_motif=num_modulate_motif/total_num

categories_revised_names=c("Create_Binding_Site", "Destroy_Binding_Site", "Modulate_Binding_Site")

pie_chart_df_prop=data.frame(matrix(NA, nrow=length(categories_revised_names), ncol=3))
colnames(pie_chart_df_prop)=c("class", "prop", "y.pos")
pie_chart_df_prop[,"prop"]=round(c(prop_create_motif, prop_destroy_motif, prop_modulate_motif), 3)
pie_chart_df_prop[,"class"]=categories_revised_names

#sort by prop
sorted_ind=order(pie_chart_df_prop[,"prop"], decreasing=TRUE)

pie_chart_df_prop=pie_chart_df_prop[sorted_ind,]

mycols <- c("green3", "gold2", "tomato1")

pie_chart_df_prop$class <- factor(pie_chart_df_prop$class, pie_chart_df_prop$class)

pie_chart_df_prop[,"y.pos"]=1-(cumsum(c(0, pie_chart_df_prop[,"prop"])) + c(pie_chart_df_prop[,"prop"] / 2, .01))[1:nrow(pie_chart_df_prop)]

p=ggplot(pie_chart_df_prop, aes(x = "", y = prop, fill = class)) +geom_bar(width = 1, stat = "identity", color = "white")+coord_polar("y", start = 0)+
	geom_label_repel(aes(y=y.pos, label = prop), size=40, show.legend = F, nudge_x = 2, nudge_y = 20)+guides(fill = guide_legend(title = ""), box.padding = unit(5, "lines"))+scale_fill_manual(values = mycols)+
	theme(legend.text=element_text(size=55, color="black"), legend.title=element_text(size=25), legend.key.size = unit(3, "cm"),  
	axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
	panel.border = element_blank(), panel.grid=element_blank(), panel.background = element_blank())

plot_name=paste(analysis_type, "binding_site", "pie_chart", sep="_")
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p, width=25, height=25)

##******************** plot out results with conservation vs without (fig. S6C) ********************

tf_alteration_analysis_without_cons_filter_cor_coef=0.25
tf_alteration_analysis_without_cons_filter_log_pvalue=-log(0.0037, base=10)
tf_alteration_analysis_with_cons_filter_cor_coef=0.37
tf_alteration_analysis_with_cons_filter_log_pvalue=-log(0.00019, base=10)

plot_df=data.frame(log_p_value=c(tf_alteration_analysis_with_cons_filter_log_pvalue, tf_alteration_analysis_without_cons_filter_log_pvalue, tf_alteration_analysis_with_cons_filter_cor_coef*10, tf_alteration_analysis_without_cons_filter_cor_coef*10), analysis_type=c("Zoonomia conservation(phyloP) filter", "No conservation filter", "Zoonomia conservation(phyloP) filter", "No conservation filter"), analysis_metric=c("-log10(p-value)", "-log10(p-value)", "correlation(r)", "correlation(r)"))

plot_df[,"analysis_metric"]=factor(as.character(plot_df[,"analysis_metric"]), c("-log10(p-value)", "correlation(r)"))
plot_df[,"analysis_type"]=factor(as.character(plot_df[,"analysis_type"]), c("Zoonomia conservation(phyloP) filter", "No conservation filter"))

mycols=c("black", "gray", "black", "gray")
p=ggplot(plot_df, aes(x = analysis_metric, y = log_p_value, fill=analysis_type))+geom_bar(width = 0.9, position = position_dodge(width = 0.9), stat = "identity")+guides(fill = guide_legend(title = "")) +scale_fill_manual(values = mycols)+xlab("")+ylab("TF motif impact vs MPRA effect")+scale_y_continuous(sec.axis = sec_axis(~ . / 10, name = ""))+theme(panel.background = element_blank(), legend.text=element_text(size=50, color="black"), legend.title=element_text(size=25), legend.key.size = unit(3, "cm"), legend.position="none", axis.title.x=element_text(size=80), axis.title.y=element_text(size=80), axis.text.x=element_text(size=50), axis.text.y=element_text(size=50), axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))

plot_name="tf_alteration_comparison_zoonomia_filter"
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p, width=15, height=20)

#******************** plot out all top hits (fig. S7) ********************#
#top_hits_df contains two columns - first column is hCONDEL_ID of interest and second column is the cell type of interest
#make your own hCONDEL list, then you can use the code to plot the MPRA results

top_hits_df=read.table("/home_path/Individual Hit Lookup/hcondel_indiv_hit_list.txt", header=F)
#top_hits_df=read.table("/home_path/Individual Hit Lookup/hcondel_indiv_hit_list.npc.txt", header=F)

colnames(top_hits_df)=c("hCONDEL_ID", "CellType")

plot_figs=TRUE
for(ind in 1:nrow(top_hits_df)){
	
	hCONDEL_ID_interest=as.character(top_hits_df[ind,1])
	cell_name=as.character(top_hits_df[ind,2])
	
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")
	lfc_skew_col=paste("log2FoldChange_Skew", cell_name, sep="_")

	chimp_padj_col=paste("padj_Chimp", cell_name, sep="_")
	human_padj_col=paste("padj_Human", cell_name, sep="_")

	lfc_chimp_col=paste("log2FoldChange_Chimp", cell_name, sep="_")
	lfc_human_col=paste("log2FoldChange_Human", cell_name, sep="_")

	lfc_se_chimp_col=paste("lfcSE_Chimp", cell_name, sep="_")
	lfc_se_human_col=paste("lfcSE_Human", cell_name, sep="_")
	
	chimp_activity_val=hcondel_full_dataset[which(hcondel_full_dataset[,"hCONDEL_ID"]==hCONDEL_ID_interest),lfc_chimp_col]
	human_activity_val=hcondel_full_dataset[which(hcondel_full_dataset[,"hCONDEL_ID"]==hCONDEL_ID_interest),lfc_human_col]

	chimp_se_val=hcondel_full_dataset[which(hcondel_full_dataset[,"hCONDEL_ID"]==hCONDEL_ID_interest),lfc_se_chimp_col]
	human_se_val=hcondel_full_dataset[which(hcondel_full_dataset[,"hCONDEL_ID"]==hCONDEL_ID_interest),lfc_se_human_col]
	
	plot_df=data.frame(value=c(chimp_activity_val, human_activity_val), sd=c(chimp_se_val, human_se_val), class=c("Chimp", "Human"))

	
	############
	
	if(plot_figs){
		mycols=c("royalblue", "limegreen")
		p=ggplot(plot_df, aes(x = class, y = value, fill=class)) +geom_bar(width = 1, stat = "identity")+geom_errorbar(aes(ymin=value-sd, ymax=value+sd), size=3)+guides(fill = guide_legend(title = "")) +scale_fill_manual(values = mycols)+xlab("")+ylab("Activity")+theme(legend.text=element_text(size=50, color="black"), legend.title=element_text(size=25), legend.key.size = unit(3, "cm"), legend.position="none", axis.title.x=element_text(size=80, margin=margin(t = 50, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=80, margin=margin(t = 0, r = 50, b = 0, l = 0)), axis.text.x=element_text(size=70), axis.text.y=element_text(size=70), panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
		plot_name=paste(ind, "MPRA_skew_barplot", sep="_")
		plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
		ggsave(filename=plot_output_file_path, plot=p, width=15, height=20)
	}
	
}


#******************** scatter plot for RNA Chimp Activity vs RNA Human Activity (Fig. 2B) ********************#

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

color_mapping=c("black", "#8DC63F", "#BE1E2D")
names(color_mapping)=c("Not Significant", "Human Gain", "Human Loss")

for(ind in 1:length(cell_names)){

	cell_name=cell_names[ind]
	
	padj_skew_cutoff=padj_val_sig_cutoff_list[[cell_name]][2]
	padj_activity_cutoff=padj_val_sig_cutoff_list[[cell_name]][1]
	
	lfc_chimp_col=paste("log2FoldChange_Chimp", cell_name, sep="_")
	lfc_human_col=paste("log2FoldChange_Human", cell_name, sep="_")
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")

	padj_human_col=paste("padj_Human", cell_name, sep="_")
	padj_chimp_col=paste("padj_Chimp", cell_name, sep="_")
	
	sig_info_skew=rep("Not Significant", nrow(hcondel_full_dataset))
	sig_info_skew[which( (hcondel_full_dataset[,lfc_human_col]<hcondel_full_dataset[,lfc_chimp_col]) & (hcondel_full_dataset[,padj_skew_col]<padj_skew_cutoff) & ( (hcondel_full_dataset[,padj_human_col]<padj_activity_cutoff) | (hcondel_full_dataset[,padj_chimp_col]<padj_activity_cutoff) ) )]="Human Loss"
	sig_info_skew[which( (hcondel_full_dataset[,lfc_human_col]>hcondel_full_dataset[,lfc_chimp_col]) & (hcondel_full_dataset[,padj_skew_col]<padj_skew_cutoff) & ( (hcondel_full_dataset[,padj_human_col]<padj_activity_cutoff) | (hcondel_full_dataset[,padj_chimp_col]<padj_activity_cutoff) ) )]="Human Gain"
	
	low_plasmid_count_indeces_human_and_chimp=low_plasmid_count_indeces_human_and_chimp_all[[cell_name]]
	
	plot_data_df=data.frame(lfc_chimp_col=hcondel_full_dataset[,lfc_chimp_col], lfc_human_col=hcondel_full_dataset[,lfc_human_col], sig_info_skew)
	
	plot_data_df=plot_data_df[-low_plasmid_count_indeces_human_and_chimp,]
	plot_data_df[,"sig_info_skew"]=factor(plot_data_df[,"sig_info_skew"], c("Not Significant", "Human Gain", "Human Loss"))
	
	plot_data_df=plot_data_df[order(plot_data_df[,"sig_info_skew"]),]
	
	plot_name=paste("log2FoldChange_Chimp_Activity", "vs", "log2FoldChange_Human_Activity", cell_name, "pval", padj_skew_cutoff, sep="_")
	p=ggplot(plot_data_df, aes(x=lfc_chimp_col, y=lfc_human_col))+xlab("log2FC Chimp Activity")+ylab("log2FC Human Activity")+
		geom_point(size=8, alpha=0.3, aes_string(colour="sig_info_skew"))+geom_vline(xintercept = 0, linetype="dashed")+geom_hline(yintercept = 0, linetype="dashed")+
		scale_color_manual(values=color_mapping)+guides(colour= guide_legend(title = "Species-Specific\nActivity"))+
		theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size=35), axis.title.y = element_text(size=35), 
		legend.title=element_text(size=30), legend.text=element_text(size=30), legend.key=element_blank(), legend.key.size = unit(1, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank())

	plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
	
	ggsave(filename=plot_output_file_path, plot=p ,width = 15, height = 15)
}

#******************** replicate correlation heatmap (fig. S5A) ********************#

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")
human_col=paste(c(rep(cell_names, each=4), rep("Plasmid", 5)), c(rep(c("R1", "R2", "R3", "R4"), 6), c("R1", "R2", "R3", "R4", "R5")), rep("Human", length(cell_names)*4+5), sep="_")
chimp_col=paste(c(rep(cell_names, each=4), rep("Plasmid", 5)), c(rep(c("R1", "R2", "R3", "R4"), 6), c("R1", "R2", "R3", "R4", "R5")), rep("Chimp", length(cell_names)*4+5), sep="_")

#combine human and chimp
hcondel_mpra_plasmid_count_table_human=hcondel_mpra_plasmid_count_table[, human_col]
colnames(hcondel_mpra_plasmid_count_table_human)=gsub("_Human", "", colnames(hcondel_mpra_plasmid_count_table_human))

hcondel_mpra_plasmid_count_table_chimp=hcondel_mpra_plasmid_count_table[, chimp_col]
colnames(hcondel_mpra_plasmid_count_table_chimp)=gsub("_Chimp", "", colnames(hcondel_mpra_plasmid_count_table_chimp))

hcondel_mpra_plasmid_count_table_human_chimp_combined=rbind(hcondel_mpra_plasmid_count_table_human, hcondel_mpra_plasmid_count_table_chimp)
hcondel_mpra_plasmid_count_table_human_chimp_combined_colnames=colnames(hcondel_mpra_plasmid_count_table_human_chimp_combined)

#change col names, remove plasmid column
keep_col=hcondel_mpra_plasmid_count_table_human_chimp_combined_colnames[which(!grepl("Plasmid", hcondel_mpra_plasmid_count_table_human_chimp_combined_colnames))]
hcondel_mpra_plasmid_count_table_iter=hcondel_mpra_plasmid_count_table_human_chimp_combined[,keep_col]

heatmap_labels=rep(cell_names, each=4)
heatmap_color_map=list("Cell Type" = c("HEK293" = "#8B1E3F", "HEPG2" = "#715470", "K562" = "#89BD9E", "GM12878" = "#F0C987", "SKNSH"="#DB4C40", "NPC"="#315C70"))
cell_name_abv=c("HEK", "HepG2", "K562", "GM12878", "SK-N-SH", "NPC")

ha = HeatmapAnnotation("Cell Type" = heatmap_labels, col = heatmap_color_map, show_legend = TRUE, show_annotation_name=FALSE, annotation_legend_param = list("Cell Type" = list(at = cell_names, labels = cell_name_abv)) )
ha_row = rowAnnotation("Cell Type" = heatmap_labels, col = heatmap_color_map, show_legend = FALSE, show_annotation_name=FALSE)
    
cor_matrix=cor(hcondel_mpra_plasmid_count_table_iter)

#remove the cell label from the correlation matrix
cell_replace_str=paste(paste(cell_names, "_", sep=""), collapse="|")
rownames(cor_matrix)=gsub(cell_replace_str, "", rownames(cor_matrix))
colnames(cor_matrix)=gsub(cell_replace_str, "", colnames(cor_matrix))

plot_name="replicate_correlation_heatmap"
plot_output_file_path=paste(c(output_file_path_with_header, plot_name, "plot.pdf"), collapse="_")

pdf(plot_output_file_path)
col_fun = colorRamp2(c(0,1) , c("white", "red"))
#include row annotation
hm=Heatmap(cor_matrix, heatmap_legend_param = list(title = "Correlation"), show_row_names = FALSE, bottom_annotation=ha, right_annotation=ha_row, column_names_rot = 90, col=col_fun, show_row_dend = FALSE, cluster_rows=FALSE, use_raster=FALSE, cluster_columns=FALSE)
draw(hm, merge_legend = TRUE)
dev.off()

#loop over and take average of correlation per cell type
all_cell_avg_cor_results=rep(NA, length(cell_names))
for(ind in 1:length(cell_names)){
	cell_name_iter=cell_names[ind]
	row_indeces=rownames(cor_matrix)[which(grepl(cell_name_iter, rownames(cor_matrix)))]
	col_indeces=colnames(cor_matrix)[which(grepl(cell_name_iter, colnames(cor_matrix)))]
	all_cell_avg_cor_results[ind]=mean(cor_matrix[row_indeces, col_indeces])
}

#take average across all cell types
mean(all_cell_avg_cor_results)

#******************** H327ac plot (fig. S5B) ********************#

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

plot_data_all=data.frame()

for(ind in 1:length(cell_names)){

	cell_name=cell_names[ind]
	lfc_human_col=paste(c("log2FoldChange_Human", cell_name), collapse="_")
	H3K27ac_col=paste("H3K27ac", cell_name, sep="_")
	
	low_plasmid_count_indeces_human=low_plasmid_count_indeces_human_all[[cell_name]]
	
	plot_data_temp=hcondel_full_dataset[-low_plasmid_count_indeces_human, c(lfc_human_col, H3K27ac_col)]
	
	colnames(plot_data_temp)=c("log2FC", "H3K27ac")
	peak_ind=which(plot_data_temp[, "H3K27ac"]==1)
	
	plot_data_temp[,"cell_name"]=cell_name
	
	plot_data_temp[,"peak_label"]="No Peak"
	plot_data_temp[peak_ind,"peak_label"]="H3K27ac Peak"

	plot_data_all=rbind(plot_data_temp, plot_data_all)
}

plot_data_all[,"cell_name"]=factor(plot_data_all[,"cell_name"], cell_names)
plot_data_all[,"peak_label"]=factor(plot_data_all[,"peak_label"], c("No Peak", "H3K27ac Peak"))

num_cell_names=length(cell_names)

y_min=-3
y_max=7
sig_info_y_start=5
plot_name="H3K27ac_log2FC_overlap"
y_axis_label="log2FC Activity"
x_axis_label=""
num_labels=2

#compare_list is the comparisons to make on the plot, compare "No Peak" vs "H3K27ac Peak" for each cell type
comparison_list=list(c("No Peak", "H3K27ac Peak"))
map_signif_level_boolean=FALSE

p=ggplot(plot_data_all, aes(x=peak_label, y=log2FC, fill=peak_label))+ylab(y_axis_label)+xlab("")+
	geom_boxplot(outlier.shape=NA)+coord_cartesian(ylim = c(y_min, y_max))+guides(fill = guide_legend(title = ""))+
	theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=30), axis.text.x=element_blank(), axis.title.y=element_text(size=30), 
	axis.text.y=element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm"), legend.key=element_blank(), 
	panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_rect(colour="black", fill = 'gray85'), 
	strip.text.x = element_text(size = 30), axis.ticks.x=element_blank())+
	facet_grid(.~cell_name, scales="free")+geom_signif(test="wilcox.test", test.args=list(alternative="less"), comparisons = comparison_list, map_signif_level=map_signif_level_boolean, tip_length=0,  size=1, textsize=8, y_position=sig_info_y_start)

plot_output_file_path=paste(output_file_path_with_header, paste(plot_name, ".pdf", sep=""), sep="_")
ggsave(filename=plot_output_file_path, plot=p, width=30, height=20 )

#********************* chimp vs human (fig. S5C) ********************#

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")
aggregated_l2fc_df=data.frame()

for(i in 1:length(cell_names)){

	cell_name=cell_names[i]
	
	low_count_indeces_chimp=low_plasmid_count_indeces_chimp_all[[cell_name]]
	low_count_indeces_human=low_plasmid_count_indeces_human_all[[cell_name]]
	
	lfc_chimp_col=paste("log2FoldChange_Chimp", cell_name, sep="_")
	lfc_human_col=paste("log2FoldChange_Human", cell_name, sep="_")

	lfc_chimp_val=hcondel_full_dataset[-low_count_indeces_chimp,lfc_chimp_col]
	l2fc_chimp_df_iter=data.frame(lfc_val=lfc_chimp_val, species=rep("Chimp", length(lfc_chimp_val)))

	lfc_human_val=hcondel_full_dataset[-low_count_indeces_human,lfc_human_col]
	l2fc_human_df_iter=data.frame(lfc_val=lfc_human_val, species=rep("Human", length(lfc_human_val)))
	
	aggregated_l2fc_df=rbind(l2fc_chimp_df_iter, l2fc_human_df_iter)
}

y_axis_label="log2FC"
x_axis_label="Species"

comparison_list=list(c("Human", "Chimp"))
map_signif_level_boolean=FALSE

sig_info_y_start=1
y_min=-2	
y_max=2

p=ggplot(aggregated_l2fc_df, aes(x=species, y=lfc_val, fill=species))+ylab(y_axis_label)+xlab(x_axis_label)+guides(fill = guide_legend(title = ""))+
	geom_boxplot(outlier.shape=NA)+geom_signif(comparisons = comparison_list, map_signif_level=map_signif_level_boolean, tip_length=0, y_position=sig_info_y_start )+coord_cartesian(ylim = c(y_min, y_max))+
	theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_text(size=25), 
	axis.ticks.x = element_blank(), legend.text=element_text(size=15), legend.key.size = unit(2, "cm"), legend.key=element_blank(), 
	panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
	
plot_name="species_vs_lfc_val_diff_categories_plot_boxplot"
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p , width=10, height=10)

#******************** allelic skew enrichments (fig. S6A) ********************#

#set analysis_type here, set to genomic features for fig. S6A, set to conservation to see enrichments with conservation
analysis_type="genomic_features"
#analysis_type="conservation"
################################################

#first get most significant skew, add to existing hcondel_full_dataset

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

padj_skew_col_names=paste("padj_Skew", cell_names, sep="_")
padj_chimp_col_names=paste("padj_Chimp", cell_names, sep="_")
padj_human_col_names=paste("padj_Human", cell_names, sep="_")

log2FoldChange_skew_col_names=paste("log2FoldChange_Skew", cell_names, sep="_")

### get other values specific to min ##

hcondel_full_dataset[,"padj_Skew_Most_Sig_Skew"]=NA
hcondel_full_dataset[,"padj_Activity_Most_Sig_Skew"]=NA
hcondel_full_dataset[,"log2FoldChange_Skew_Most_Sig_Skew"]=NA

hcondel_full_dataset[,"Most_Sig_Skew_abs_predicted_JASPARv2020_tf_hg38_panTro4_diff"]=NA
hcondel_full_dataset[,"Most_Sig_Skew_H3K27ac"]=NA

H3K27ac_cols=paste("H3K27ac", cell_names, sep="_")

tf_diff_cols=paste(cell_names, "predicted_JASPARv2020_tf_hg38_panTro4_diff", sep="_")

for(ind in 1:nrow(hcondel_full_dataset)){

	### skew ###
	min_skew_ind=which.min(hcondel_full_dataset[ind,padj_skew_col_names])
	
	if(length(min_skew_ind)!=0){
		
		padj_chimp_iter=hcondel_full_dataset[ind,padj_chimp_col_names[min_skew_ind]]
		padj_human_iter=hcondel_full_dataset[ind,padj_human_col_names[min_skew_ind]]
		padj_skew_iter=hcondel_full_dataset[ind,padj_skew_col_names[min_skew_ind]]
	
		log2FoldChange_skew_iter=hcondel_full_dataset[ind,log2FoldChange_skew_col_names[min_skew_ind]]
	
		padj_activity_min_matched=min(padj_chimp_iter, padj_human_iter)
	
		hcondel_full_dataset[ind,"padj_Skew_Most_Sig_Skew"]=padj_skew_iter
		hcondel_full_dataset[ind,"padj_Activity_Most_Sig_Skew"]=padj_activity_min_matched
		hcondel_full_dataset[ind,"log2FoldChange_Skew_Most_Sig_Skew"]=log2FoldChange_skew_iter
	
		#get the H3K27ac signal and JASPAR abs TF binding difference signal
		hcondel_full_dataset[ind,"Most_Sig_Skew_abs_predicted_JASPARv2020_tf_hg38_panTro4_diff"]=abs(hcondel_full_dataset[ind,tf_diff_cols[min_skew_ind]])
		hcondel_full_dataset[ind,"Most_Sig_Skew_H3K27ac"]=hcondel_full_dataset[ind,H3K27ac_cols[min_skew_ind]]
	}
	
}

hcondel_full_dataset[,"log_padj_Activity_Most_Sig_Skew"]=-log(hcondel_full_dataset[,"padj_Activity_Most_Sig_Skew"], base=10)
hcondel_full_dataset[which(hcondel_full_dataset[,"log_padj_Activity_Most_Sig_Skew"]==Inf), "log_padj_Activity_Most_Sig_Skew"]=10^300

hcondel_full_dataset[,"del_len"]=hcondel_full_dataset[,"panTro4_del_end_pos"]-hcondel_full_dataset[,"panTro4_del_start_pos"]

#set positions without any TF binding/alteration scores to 0
hcondel_full_dataset[which(is.na(hcondel_full_dataset[,"Most_Sig_Skew_abs_predicted_JASPARv2020_tf_hg38_panTro4_diff"])),"Most_Sig_Skew_abs_predicted_JASPARv2020_tf_hg38_panTro4_diff"]=0

#now run allelic skew enrichments

gen_allelic_skew_enrichment_stats=function(data_table, var_of_interest_info_table, threshold_info_table){

	#need to loop over type_names to determine whether or not to make aggregated plot
	enrichment_output_pval=rep(NA, length(var_of_interest_info_table[,1]))
	enrichment_output_coef=rep(NA, length(var_of_interest_info_table[,1]))
	
	names(enrichment_output_pval)=var_of_interest_info_table[,1]
	names(enrichment_output_coef)=var_of_interest_info_table[,1]
		
	for(i in 1:nrow(var_of_interest_info_table)){
	
		var_interest=as.character(var_of_interest_info_table[i,1])
		var_control=as.character(var_of_interest_info_table[i,3])
		test_direction=as.character(var_of_interest_info_table[i,4])
		
		item_label=rep(0, nrow(data_table))

		#mark any indeces that passes the cutoff values given
		threshold_info_table_names=names(threshold_info_table)
		
		threshold_info_boolean_table=data.frame(matrix(FALSE, nrow=nrow(data_table), ncol=length(threshold_info_table_names)))
		
		for (j in 1:length(threshold_info_table_names)){
		
			threshold_name=threshold_info_table_names[j]
			threshold_value=threshold_info_table[[threshold_name]]
			
			pass_ind=which(data_table[,threshold_name]<threshold_value)
			threshold_info_boolean_table[pass_ind,j]=TRUE
		}
		
		all_pass_ind=which(apply(threshold_info_boolean_table, 1, all))
		item_label[all_pass_ind]=1

		#combine label with values 
		lm_data=data.frame(item_label, value=data_table[,var_interest], control_value=data_table[,var_control])
		
		#run regression here
		lm_output=lm(item_label ~ value + control_value, lm_data)
		
		p_val=summary(lm_output)$coefficients[2,4]
		coef=summary(lm_output)$coefficients[2,1]
		
		#store the log pvalues for plotting later
		enrichment_output_pval[i]=-log(p_val, base=10)
		enrichment_output_coef[i]=coef
	}			

	return(list(enrichment_output_pval, enrichment_output_coef))
}


gen_comparison_barplots=function(gen_compare_dist_plots_result, var_of_interest_info_table, save_file_path_header, fdr_cutoff){

	enrichment_output_pval=gen_compare_dist_plots_result[[1]]
	enrichment_output_coef=gen_compare_dist_plots_result[[2]]

	var_names=names(enrichment_output_pval)
	FDR_info=rep("",length(enrichment_output_pval))
	p_adj_vals=p.adjust(10^(-enrichment_output_pval), method="BH")
	FDR_pass_ind=which(p_adj_vals<fdr_cutoff)
	FDR_info[FDR_pass_ind]="*"

	plot_data=data.frame(FDR_info=FDR_info, variable=var_names, log_pval=enrichment_output_pval, coef=enrichment_output_coef)
	
	########### change long column names #######
	
	plot_data[,"variable"]=as.character(plot_data[,"variable"])
	for (i in 1:nrow(var_of_interest_info_table)){
		current_name_info=var_of_interest_info_table[i,1]
		new_name_info=var_of_interest_info_table[i,2]
		
		plot_data[which(plot_data[,"variable"]==current_name_info),"variable"]=new_name_info
	}
	
	########### plot pvals ####################################
	
	plot_data=plot_data[order(plot_data[,"log_pval"]), ]
	plot_data[,"variable"]=factor(plot_data[,"variable"], plot_data[,"variable"])
	
	plot_data[,"coef_directionality"]="Negative"
	pos_ind=which(plot_data[,"coef"]>0)
	plot_data[pos_ind,"coef_directionality"]="Positive"
	plot_data[,"coef_directionality"]=factor(plot_data[,"coef_directionality"], c("Negative", "Positive"))
	
	#fill with coef directionality
	#plot stars for significance
	mycols=c("royalblue", "limegreen")
	p=ggplot(plot_data, aes(x = variable, y=log_pval, fill=coef_directionality)) +geom_bar(stat="identity")+geom_text(aes(label=FDR_info), size=20, hjust=-0.25)+
		xlab("")+ylab("Coefficient log(p-value)")+scale_fill_manual("Coefficient\nDirectionality", values=mycols)+ 
		theme(legend.text=element_text(size=30, color="black"), legend.title=element_text(size=30, color="black"), legend.key.size = unit(2, "cm"), plot.title = element_text(hjust = 0.5), 
		axis.text.x=element_text(size=30), axis.text.y=element_text(size=30), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), 
		panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_flip()
	
	plotName="lm_pval_barplot"
	save_file_path=paste(c(save_file_path_header,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
	
	ggsave(filename=save_file_path, plot=p, width = 20, height = 20, limitsize=FALSE)

	######### plot coefficients #############
	
	plot_data=plot_data[order(plot_data[,"coef"]), ]
	plot_data[,"variable"]=factor(plot_data[,"variable"], plot_data[,"variable"])
	
	p=ggplot(plot_data, aes(x = variable, y=coef, fill=FDR_info)) +geom_bar(stat="identity")+xlab("")+ylab("Coefficient Value")+
		theme(legend.text=element_text(size=20, color="black"), legend.title=element_text(size=20), legend.key.size = unit(2, "cm"), plot.title = element_text(hjust = 0.5), 
		axis.text.x=element_text(size=30), axis.text.y=element_text(size=30), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), 
		panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_flip()

	plotName="lm_coef_barplot"
	save_file_path=paste(c(save_file_path_header,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
	
	ggsave(filename=save_file_path, plot=p, width = 20, height = 20, limitsize=FALSE)

}


#filter out oligos with too many mismatches/unaligned
num_bp_unaligned_cutoff=quantile(hcondel_full_dataset[,"num_bp_unaligned"], probs = seq(0, 1, 0.1))[10]
num_bp_mismatches_cutoff=quantile(hcondel_full_dataset[,"num_bp_mismatches"], probs = seq(0, 1, 0.1))[9]

plasmid_cutoff_val=60

filter_ind=which( (hcondel_full_dataset[,"plasmid_avg_chimp"]<plasmid_cutoff_val) | (hcondel_full_dataset[,"plasmid_avg_human"]<plasmid_cutoff_val) | (hcondel_full_dataset[,"num_bp_unaligned"]>=num_bp_unaligned_cutoff) | (hcondel_full_dataset[,"num_bp_mismatches"]>=num_bp_mismatches_cutoff))
hcondel_full_dataset_filtered=hcondel_full_dataset[-filter_ind, ]

var_control="log_padj_Activity_Most_Sig_Skew"
#create threshold_info_table, value is cutoff (i.e. padj value cutoff)
threshold_info_table=list()
threshold_info_table[["padj_Skew_Most_Sig_Skew"]]=0.2
#threshold_info_table[["padj_Activity_Most_Sig_Skew"]]=0.1

if(analysis_type=="genomic_features"){

	#first column is variable name, second column is name for plot, third column is test direction
	var_of_interest_info_table=var_of_interest_info_table=data.frame(matrix(nrow=3, ncol=4))
	var_of_interest_info_table[1,]=c("del_phyloP_max_score", "Zoonomia phyloP Max", var_control, "greater")
	var_of_interest_info_table[2,]=c("ENCODE_ccre", "Screen CCRE", var_control, "greater")
	var_of_interest_info_table[3,]=c("Most_Sig_Skew_abs_predicted_JASPARv2020_tf_hg38_panTro4_diff", "TF Binding Difference", var_control, "greater")
	
} else if(analysis_type=="conservation"){

	#first column is variable name, second column is name for plot, third column is test direction
	var_of_interest_info_table=var_of_interest_info_table=data.frame(matrix(nrow=2, ncol=4))
	var_of_interest_info_table[1,]=c("cons_phastCons_score", "Conserved Block PhastCons Score", var_control, "greater")
	var_of_interest_info_table[2,]=c("del_phyloP_max_score", "Zoonomia phyloP Max", var_control, "greater")
} else{
	print("Invalid Analysis Type!")
}

##################################################################

gen_allelic_skew_enrichment_stats_results=gen_allelic_skew_enrichment_stats(hcondel_full_dataset_filtered, var_of_interest_info_table, threshold_info_table)

# take previous output and create bar plots

save_file_path_header=paste(output_file_path_with_header, analysis_type, sep="_")
fdr_cutoff=0.05

gen_comparison_barplots(gen_allelic_skew_enrichment_stats_results, var_of_interest_info_table, save_file_path_header, fdr_cutoff)


print("pval")
10^(-gen_allelic_skew_enrichment_stats_results[[1]])

print("coefficient")
gen_allelic_skew_enrichment_stats_results[[2]]

#******************** LOXL2 HCR-FlowFish results (fig. S8B) ********************#

#obtain count data from crispresso analyses
rep_1_vals=c(142489, 120718, 51534, 131736)
rep_2_vals=c(221351, 94356, 55701, 47385)

count_table=matrix(NA, nrow=2, ncol=4)
count_table[1,]=rep_1_vals
count_table[2,]=rep_2_vals

iter_analysis_table=data.frame(matrix(nrow=2, ncol=9))
colnames(iter_analysis_table)=c("id", "odds_ratio", "p_val", "odds_ratio_lb", "odds_ratio_ub", "human_high", "human_low", "chimp_high", "chimp_low")

#run fisher test, get stats
num_replicates=nrow(count_table)

for(ind in 1:num_replicates){

	count_values=count_table[ind,]
	fisher_table=matrix(count_values, nrow=2, ncol=2)
	
	fisher_test_results=fisher.test(fisher_table)
	p_val=fisher_test_results$p.value
	odds_est=fisher_test_results$estimate
	odds_est_lb=fisher_test_results$conf.int[[1]]
	odds_est_ub=fisher_test_results$conf.int[[2]]

	iter_analysis_table[ind,2:5]=c(odds_est, p_val, odds_est_lb, odds_est_ub)
	iter_analysis_table[ind,6:9]=count_values
}

iter_analysis_table[, "id"]=paste("Rep", seq(1, num_replicates), sep=" ")
iter_analysis_table[, "odds_ratio"]=log(iter_analysis_table[, "odds_ratio"], base=2)
iter_analysis_table[, "odds_ratio_lb"]=log(iter_analysis_table[, "odds_ratio_lb"], base=2)
iter_analysis_table[, "odds_ratio_ub"]=log(iter_analysis_table[, "odds_ratio_ub"], base=2)

ymin_plot=-0.3
ymax_plot=2

p=ggplot(iter_analysis_table, aes(x = id, y = odds_ratio)) +geom_bar(width = 1, stat = "identity")+geom_errorbar(aes(ymin=odds_ratio_lb, ymax=odds_ratio_ub), size=3)+
	geom_hline(aes(yintercept=0), size=0.7, color='red', linetype="dashed")+guides(fill = guide_legend(title = ""))+xlab("")+ylab("log2FC ratio")+
	theme(panel.background = element_blank(), legend.text=element_text(size=50, color="black"), legend.title=element_text(size=25), legend.key.size = unit(3, "cm"), legend.position="none", 
	axis.title.x=element_text(size=80), axis.title.y=element_text(size=80), axis.text.x=element_text(size=70), axis.text.y=element_text(size=70), 
	axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))+ylim(c(ymin_plot, ymax_plot))

plot_name="LOXL2_HCR_results_barplot"
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p, width=25, height=20)

#******************** MISCELLANEOUS plots/statistics for paper ********************#

#******************** plotting count distributions for significant hits ********************#
library(fancycut)

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

plasmid_cutoff_map=list()
plasmid_cutoff_map[["HEK293"]]=20
plasmid_cutoff_map[["HEPG2"]]=20
plasmid_cutoff_map[["K562"]]=20
plasmid_cutoff_map[["GM12878"]]=20
plasmid_cutoff_map[["SKNSH"]]=20
plasmid_cutoff_map[["NPC"]]=60

hcondel_full_dataset[,"plasmid_avg_human_chimp"]=(hcondel_full_dataset[,"plasmid_avg_human"]+hcondel_full_dataset[,"plasmid_avg_chimp"])/2

num_count_bins=10

cell_plasmid_count_info_table_keep_all=data.frame()

for(i in 1:length(cell_names)){
	cell_name=cell_names[i]
	
	plasmid_cutoff_val=plasmid_cutoff_map[[cell_name]]
	
	#make quantile bins, then filter out low plasmid counts
	
	plasmid_counts=hcondel_full_dataset[, "plasmid_avg_human_chimp"]

	#get quantile labels
	quantile_cutoffs=quantile(plasmid_counts, seq(0, 1, 1/num_count_bins))

	interval_vec=rep(NA, num_count_bins)
	plasmids_bin_categories_names=rep(NA, num_count_bins) 
	for(ind in 2:length(quantile_cutoffs)){
		interval=paste('(', quantile_cutoffs[ind-1], ',', quantile_cutoffs[ind], ']',  sep="")
		interval_vec[ind-1]=interval
		plasmids_bin_categories_names[ind-1]=paste(names(quantile_cutoffs)[ind-1], names(quantile_cutoffs)[ind], sep=" - ")
	}

	plasmids_bin_categories=paste(interval_vec, sep=":")

	plasmid_count_bin_labels=wafflecut(plasmid_counts, plasmids_bin_categories, plasmids_bin_categories_names, na.bucket="missing")

	cell_plasmid_count_info_table=data.frame(hCONDEL_ID=as.character(hcondel_full_dataset[, "hCONDEL_ID"]), plasmid_avg_human=hcondel_full_dataset[, "plasmid_avg_human"], plasmid_avg_chimp=hcondel_full_dataset[, "plasmid_avg_chimp"], plasmid_avg_human_chimp=hcondel_full_dataset[, "plasmid_avg_human_chimp"], quantile_bin=plasmid_count_bin_labels)

	cell_plasmid_count_info_table[,"quantile_bin"]=factor(as.character(cell_plasmid_count_info_table[,"quantile_bin"]), plasmids_bin_categories_names)
	
	cell_plasmid_count_info_table[,"type"]=rep(cell_name, length(plasmid_counts))
	
	#filter out low plasmid counts
	cell_plasmid_count_info_table=cell_plasmid_count_info_table[-which( (cell_plasmid_count_info_table[, "plasmid_avg_human"]<plasmid_cutoff_val) | (cell_plasmid_count_info_table[, "plasmid_avg_chimp"]<plasmid_cutoff_val) ),]
	
	padj_skew_cutoff=padj_val_sig_cutoff_list[[cell_name]][2]
	padj_activity_cutoff=padj_val_sig_cutoff_list[[cell_name]][1]
	
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")

	padj_human_col=paste("padj_Human", cell_name, sep="_")
	padj_chimp_col=paste("padj_Chimp", cell_name, sep="_")
	
	#classify hits into each quantile bin
	sig_seq_ids=as.character(hcondel_full_dataset[which( (hcondel_full_dataset[,padj_skew_col]<padj_skew_cutoff) & ( (hcondel_full_dataset[,padj_human_col]<padj_activity_cutoff) | (hcondel_full_dataset[,padj_chimp_col]<padj_activity_cutoff) ) ), "hCONDEL_ID"])
	
	sig_ind=which(as.character(cell_plasmid_count_info_table[, "hCONDEL_ID"]) %in% sig_seq_ids)
	
	cell_plasmid_count_info_table_keep=cell_plasmid_count_info_table[sig_ind,]
	
	print(cell_name)
	print("proportion of skew hits in lowest 10% plasmid count bin:")
	print(length(which(as.character(cell_plasmid_count_info_table_keep[,"quantile_bin"])==plasmids_bin_categories_names[1]))/nrow(cell_plasmid_count_info_table_keep))
	
	#num sig hits
	#print(nrow(cell_plasmid_count_info_table_keep))
	cell_plasmid_count_info_table_keep_all=rbind(cell_plasmid_count_info_table_keep_all, cell_plasmid_count_info_table_keep)	
}

#plot proportions

plot_name="sig_skew_plasmid_count_quantile_proportion"
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_plasmid_count_info_table_keep_all, aes_string(x="type", fill="quantile_bin"))+geom_bar(position="fill")+xlab("Cell Type")+ylab("Proportion of hCONDEL skew hits in quantile bin")+scale_fill_discrete(name="Plasmid count bin percent")+theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=25), axis.title.y=element_text(size=40), axis.text.y=element_text(size=25), legend.title = element_text(size=30), legend.text=element_text(size=25), legend.key.size = unit(1.5, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(filename=plot_output_file_path, plot=p ,width = 20, height = 20)

#******************** look at plasmid count cutoffs on skew sd  ********************#

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

plasmid_cutoff_map=list()
plasmid_cutoff_map[["HEK293"]]=20
plasmid_cutoff_map[["HEPG2"]]=20
plasmid_cutoff_map[["K562"]]=20
plasmid_cutoff_map[["GM12878"]]=20
plasmid_cutoff_map[["SKNSH"]]=20
plasmid_cutoff_map[["NPC"]]=60

hcondel_full_dataset[,"plasmid_avg_human_chimp"]=(hcondel_full_dataset[,"plasmid_avg_human"]+hcondel_full_dataset[,"plasmid_avg_chimp"])/2

for(j in 1:length(cell_names)){
	cell_name=cell_names[j]

	plasmid_cutoff_val=plasmid_cutoff_map[[cell_name]]
	
	padj_skew_cutoff=padj_val_sig_cutoff_list[[cell_name]][2]
	padj_activity_cutoff=padj_val_sig_cutoff_list[[cell_name]][1]
	
	lfcSE_skew_col=paste("lfcSE_Skew", cell_name, sep="_")
	lfc_skew_col=paste("log2FoldChange_Skew", cell_name, sep="_")
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")
	
	padj_human_col=paste("padj_Human", cell_name, sep="_")
	padj_chimp_col=paste("padj_Chimp", cell_name, sep="_")
	
	hcondel_full_dataset_plasmid_count_iter=hcondel_full_dataset[,c(lfc_skew_col, lfcSE_skew_col, padj_skew_col, padj_human_col, padj_chimp_col, "plasmid_avg_human", "plasmid_avg_chimp", "plasmid_avg_human_chimp")]
	
	plasmid_high_count_ind=which( (hcondel_full_dataset_plasmid_count_iter[,"plasmid_avg_human"]>=plasmid_cutoff_val) & (hcondel_full_dataset_plasmid_count_iter[,"plasmid_avg_chimp"]>=plasmid_cutoff_val))
	
	#how many hCONDELs retained from plasmid filter
	#print(length(plasmid_high_count_ind))
	
	low_count_label=paste("Chimp or Human plasmid count \n < than ", plasmid_cutoff_val, sep="")
	high_count_label=paste("Chimp and Human plasmid count \n >= than ", plasmid_cutoff_val, sep="")
	
	hcondel_full_dataset_plasmid_count_iter[, "plasmid_info"]=low_count_label
	hcondel_full_dataset_plasmid_count_iter[plasmid_high_count_ind, "plasmid_info"]=high_count_label
	
	hcondel_full_dataset_plasmid_count_iter_sig=hcondel_full_dataset_plasmid_count_iter[which( (hcondel_full_dataset_plasmid_count_iter[,padj_skew_col]<padj_skew_cutoff) & ( (hcondel_full_dataset_plasmid_count_iter[,padj_human_col]<padj_activity_cutoff) | (hcondel_full_dataset_plasmid_count_iter[,padj_chimp_col]<padj_activity_cutoff) ) ), ]

	mycols = list()
	mycols[low_count_label]="red"
	mycols[high_count_label]="forestgreen"
	
	#print(table(hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_info"])[1]/nrow(hcondel_full_dataset_plasmid_count_iter_sig))
	
	#low_quantile_count=quantile(hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_avg_human_chimp"], 0.1)
	#prop_high_skew_in_low_quantile=length(which( (hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_avg_human_chimp"]<low_quantile_count) &  (abs(hcondel_full_dataset_plasmid_count_iter_sig[,lfc_skew_col])>2) ))/length(which(hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_avg_human_chimp"]<low_quantile_count))
	#print(prop_high_skew_in_low_quantile)

	low_quantile_count=quantile(hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_avg_human_chimp"], 0.05)
	prop_high_skew_in_low_quantile=length(which( (hcondel_full_dataset_plasmid_count_iter_sig[,"plasmid_avg_human_chimp"]<low_quantile_count) &  (abs(hcondel_full_dataset_plasmid_count_iter_sig[,lfc_skew_col])>1) ))/length(which(abs(hcondel_full_dataset_plasmid_count_iter_sig[,lfc_skew_col])>1))
	print(prop_high_skew_in_low_quantile)
		
	#plot points for all cell types
	plot_name=paste(cell_name, "lfc_skew", "padj_cutoff", padj_skew_cutoff, "plasmid_cutoff", plasmid_cutoff_val, sep="_")
	p=ggplot(hcondel_full_dataset_plasmid_count_iter_sig, aes_string(x="plasmid_avg_human_chimp", y=lfc_skew_col))+xlab("Average plasmid count")+ylab(paste("l2fc skew"))+geom_point(size=8, alpha=0.3, aes_string(colour="plasmid_info"))+scale_color_manual(values=mycols)+guides(colour= guide_legend(title = ""))+theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size=35), axis.title.y = element_text(size=35), legend.key=element_blank(), legend.text=element_text(size=25), legend.key.size = unit(2, "cm"), legend.position = c(0.87, 0.25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
	plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
	ggsave(filename=plot_output_file_path, plot=p ,width = 20, height = 20)	
	
	#plot points for all cell types
	plot_name=paste(cell_name, "lfcSE_skew", "padj_cutoff", padj_skew_cutoff, "plasmid_cutoff", plasmid_cutoff_val, sep="_")
	p=ggplot(hcondel_full_dataset_plasmid_count_iter_sig, aes_string(x="plasmid_avg_human_chimp", y=lfcSE_skew_col))+xlab("Average plasmid count")+ylab(paste("lfcSE skew"))+geom_point(size=8, alpha=0.3, aes_string(colour="plasmid_info"))+scale_color_manual(values=mycols)+guides(colour= guide_legend(title = ""))+theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size=35), axis.title.y = element_text(size=35), legend.key=element_blank(), legend.text=element_text(size=25), legend.key.size = unit(2, "cm"), legend.position = c(0.87, 0.25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
	plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
	ggsave(filename=plot_output_file_path, plot=p ,width = 20, height = 20)	

}

#******************** skew sd for each cell type for significant hits ********************#
#plot skew sd all or significant skew only
#analysis_type="all_hcondels"
analysis_type="significant_skew"

cell_names=c("HEK293", "HEPG2", "K562", "GM12878", "SKNSH", "NPC")

padj_val_sig_cutoff_list=list()
padj_val_sig_cutoff_list[["HEK293"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["HEPG2"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["K562"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["GM12878"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["SKNSH"]]=c(0.1, 0.05)
padj_val_sig_cutoff_list[["NPC"]]=c(0.1, 0.05)

plot_data=data.frame()

for(j in 1:length(cell_names)){
	cell_name=cell_names[j]
	
	padj_skew_cutoff=padj_val_sig_cutoff_list[[paste(cell_name, sep="_")]][2]

	lfcSE_skew_col=paste("lfcSE_Skew", cell_name, sep="_")
	lfc_skew_col=paste("log2FoldChange_Skew", cell_name, sep="_")
	padj_skew_col=paste("padj_Skew", cell_name, sep="_")
	
	padj_human_col=paste("padj_Human", cell_name, sep="_")
	padj_chimp_col=paste("padj_Chimp", cell_name, sep="_")
			
	if(analysis_type=="all_hcondels"){
		hcondel_full_dataset_filtered=hcondel_full_dataset[, c(lfc_skew_col, lfcSE_skew_col, padj_skew_col)]
	} else if(analysis_type=="significant_skew"){
		low_plasmid_count_indeces_human_and_chimp=low_plasmid_count_indeces_human_and_chimp_all[[cell_name]]
		hcondel_full_dataset_filtered=hcondel_full_dataset[-low_plasmid_count_indeces_human_and_chimp, c(lfc_skew_col, lfcSE_skew_col, padj_skew_col)]
		hcondel_full_dataset_filtered=hcondel_full_dataset_filtered[which( (hcondel_full_dataset[,padj_skew_col]<padj_skew_cutoff) & ( (hcondel_full_dataset[,padj_human_col]<padj_activity_cutoff) | (hcondel_full_dataset[,padj_chimp_col]<padj_activity_cutoff) ) ),]
	} else{
		print("Invalid Analysis Type!")
		break
	}
	
	colnames(hcondel_full_dataset_filtered)=gsub(paste("", cell_name, sep="_"), "", colnames(hcondel_full_dataset_filtered))
	hcondel_full_dataset_filtered[,"type"]=cell_name
	
	plot_data=rbind(plot_data, hcondel_full_dataset_filtered)	
}

#mycols = list("HEK293FT" = "#8B1E3F", "HEPG2" = "#715470", "K562" = "#89BD9E", "GM12878" = "#F0C987", "SKNSH"="#DB4C40", "NPC"="#315C70")

p=ggplot(plot_data, aes(x=type, y=lfcSE_Skew))+ylab("log2FC skew standard error")+xlab("Cell type")+geom_boxplot()+theme(axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#plot points for all cell types
plot_name=paste("lfcSE_skew_all_cell_types", analysis_type, sep="_")
plot_output_file_path=paste(c(output_file_path_with_header,paste(c("_",plot_name),collapse=""), ".pdf"), collapse="")
ggsave(filename=plot_output_file_path, plot=p ,width = 20, height = 20)

