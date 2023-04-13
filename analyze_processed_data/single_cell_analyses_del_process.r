library(dplyr)
library(Seurat)
library(patchwork)
library(goseq)
library(ggplot2)
library(ggsignif)
library(EnhancedVolcano)

setwd("/home_path/10X")
out_file_path="/home_path/10X/Analysis_Output"
out_prefix="del_1a_2_13_21"

out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")

#******************** preprocess single cell data ********************#

#load data independently, then combine
rep1.data <- Read10X(data.dir = "/home_path/10X/del_1a_2_13_21/rep1_rna_out/outs/filtered_feature_bc_matrix/")
rep1_info <- CreateSeuratObject(counts = rep1.data, project = "rep1")
 
rep2.data <- Read10X(data.dir = "/home_path/10X/del_1a_2_13_21/rep2_rna_out/outs/filtered_feature_bc_matrix/")
rep2_info <- CreateSeuratObject(counts = rep2.data, project = "rep2")
 
cell_info <- merge(x = rep1_info, y = rep2_info, add.cell.ids=c("rep1", "rep2"), project = "del_1a_2_13_21_combined")
                   
cell_info

#get pct mt genes
cell_info[["percent.mt"]] <- PercentageFeatureSet(cell_info, pattern = "^MT-")

#visualize QC
VlnPlot(cell_info, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(cell_info, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cell_info, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#cell_info <- subset(cell_info, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#remove dead cells
cell_info <- subset(cell_info, subset =  percent.mt < 19)

#data normalization, 10000 is default
cell_info <- NormalizeData(cell_info, normalization.method = "LogNormalize", scale.factor = 10000)

#find most variable genes

cell_info <- FindVariableFeatures(cell_info, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cell_info), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cell_info)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling data
all.genes <- rownames(cell_info)
cell_info <- ScaleData(cell_info, features = all.genes)

#linear reduction
cell_info <- RunPCA(cell_info, features = VariableFeatures(object = cell_info))

# Examine and visualize PCA results a few different ways
print(cell_info[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(cell_info, dims = 1:2, reduction = "pca")

DimPlot(cell_info, reduction = "pca")
DimHeatmap(cell_info, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#cell_info <- JackStraw(cell_info, num.replicate = 100)
#cell_info <- ScoreJackStraw(cell_info, dims = 1:20)

#JackStrawPlot(cell_info, dims = 1:15)

#elbow plot
ElbowPlot(cell_info)

#clustering cells
cell_info <- FindNeighbors(cell_info, dims = 1:10)
cell_info <- FindClusters(cell_info, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(cell_info), 5)

#run UMAP

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
cell_info <- RunUMAP(cell_info, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(cell_info, reduction = "umap")

############## partition cell types into more categories ###################

### rep ###
cell_info@meta.data$rep_info=cell_info@meta.data$orig.ident

### N vs S ###
#determined by looking at seurat clusters beforehand (see below for plotting)
cell_names=rownames(cell_info@meta.data)

clust_1_ind=which(cell_info@meta.data$seurat_clusters%in%as.character(c(1,2,3,4,5,7)))
cell_names_clust_1=cell_names[clust_1_ind]

clust_2_ind=which(cell_info@meta.data$seurat_clusters%in%as.character(c(0,6,8,9,10,11,12)))
cell_names_clust_2=cell_names[clust_2_ind]

cell_info@meta.data$cell_type=rep(NA, nrow(cell_info@meta.data))
cell_info@meta.data$cell_type[clust_1_ind]="S"

cell_info@meta.data$cell_type[clust_2_ind]="N"


#### phase #######

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cell_info <- CellCycleScoring(cell_info, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(cell_info, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
cell_info <- RunPCA(cell_info, features = c(s.genes, g2m.genes))

#for changing labels
#https://satijalab.org/seurat/v3.0/interaction_vignette.html


#### hdr ids ###

hdr_info_file_rep1="/home_path/10X/del_1a_2_13_21/rep1_targeted_loxl2_rna_parallel.summTable.concat.umi_collapsed.id.out.txt"
hdr_info_df_rep1=read.table(hdr_info_file_rep1, header=T)

#modify id
hdr_info_df_rep1[,1]=paste("rep1_", hdr_info_df_rep1[,1] , "-1", sep="")

hdr_info_file_rep2="/home_path/10X/del_1a_2_13_21/rep2_targeted_loxl2_rna_parallel.summTable.concat.umi_collapsed.id.out.txt"
hdr_info_df_rep2=read.table(hdr_info_file_rep2, header=T)

#modify id
hdr_info_df_rep2[,1]=paste("rep2_", hdr_info_df_rep2[,1] , "-2", sep="")

hdr_info_df=rbind(hdr_info_df_rep1, hdr_info_df_rep2)

matched_ind=match(rownames(cell_info@meta.data), hdr_info_df[,"BC"])
#length(which(!is.na(matched_ind)))

cell_info@meta.data$hdr_info=hdr_info_df[matched_ind, "Type"]

#**************** some diagnostic plots ****************#

#### diagnostics ###

plotName="phase"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
DimPlot(cell_info, pt.size=1, group.by="Phase", shuffle=TRUE)
dev.off()

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score")

plotName="metrics"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
FeaturePlot(cell_info, reduction = "umap", features = metrics, pt.size = 0.4, sort.cell = TRUE,min.cutoff = 'q10')
dev.off()

plotName="rep"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
DimPlot(cell_info, group.by="rep_info")
dev.off()


#look at seurat clusters (helps see which subclusters belong to S, which one belongs to N)
plotName="seurat_clusters"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
DimPlot(cell_info, pt.size=1, group.by="seurat_clusters", shuffle=TRUE)
dev.off()

#CD44 to distinguish S vs N
#CD44 is characteristically expressed in S-cells
#10.1593/neo.04310
#10.1007/s13402-011-0022-z

marker_genes=c("CD44")

plotName="hdr_genes_CD44"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
FeaturePlot(cell_info, reduction = "umap", features = marker_genes, sort.cell = TRUE, min.cutoff = 'q10')
dev.off()

#color by specific marker genes as based on fig 1H in https://doi.org/10.1038/ng.3921, also helps to distinguish S vs N
marker_genes=c("HAND2", "PHOX2A", "PHOX2B", "GATA2", "GATA3", "HAND1", "KLF7", "ISL1")

plotName="marker_1_genes"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
FeaturePlot(cell_info, reduction = "umap", features = marker_genes, sort.cell = TRUE, min.cutoff = 'q10')
dev.off()

marker_genes=c("NR3C1", "BHLHE41", "MAFF", "GLIS3", "IRF1", "IRF2", "IRF3", "FLI1", "MEF2D", "PRRX1", "RUNX1", "RUNX2", "TBX18", "FOSL1", "FOSL2")

plotName="marker_2_genes"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

pdf(plotOutFileName, width=15, height=10)
par(cex.axis=1, cex.lab=1, cex.main=1, cex.sub=1)
FeaturePlot(cell_info, reduction = "umap", features = marker_genes, sort.cell = TRUE, min.cutoff = 'q10')
dev.off()

#### checking LOXL2 expression in entire UMAP #####
marker_genes=c("LOXL2")

plotName="hdr_genes_LOXL2"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=FeaturePlot(cell_info, reduction = "umap", features = marker_genes, pt.size = 4, sort.cell = TRUE, cols=c("lightgrey", "red"), min.cutoff = 'q10')+theme(plot.title = element_text(size=40), axis.text.x = element_text(size=40), axis.text.y = element_text(size=40), axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.text = element_text(size=40), legend.key.size = unit(2, 'cm'))
ggsave(filename=plotOutFileName, plot=p, width=25, height=20 )


#loxl2 plot showing low expression in N, high in S
plotName="loxl2_vln_plot_S_vs_N"
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")

cell_info <- SetIdent(cell_info, value = "hdr_info")
p=VlnPlot(cell_info, features=c("LOXL2"), idents=c("WT", "HDR"), group.by="cell_type", cols=c("blue", "orange"), pt.size=3)+theme(plot.title = element_text(size=40), axis.text.x = element_text(size=40), axis.text.y = element_text(size=40), axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.text = element_text(size=40), legend.position = "none")
ggsave(filename=plotOutFileName, plot=p, width=25, height=20 )


#**************** Fig. 4D,E plots (UMAP, LOXL2 violin expression) ****************#
#only use S-type cells as those have LOXL2 expressed

cell_info_parsed_df=as.data.frame(cell_info[["umap"]]@cell.embeddings)

expression_data=GetAssayData(object = cell_info, slot = "data")["LOXL2",]

cell_info_parsed_df[,"exp"]=as.numeric(expression_data)

cell_info_parsed_df[,"cell_type"]=cell_info@meta.data$cell_type
cell_info_parsed_df[,"hdr_info"]=cell_info@meta.data$hdr_info

#sanity checks
all(names(expression_data)==rownames(cell_info_parsed_df))
all(names(expression_data)==rownames(cell_info@meta.data))

cell_info_parsed_df_s_wt_hdr_only=cell_info_parsed_df[which( (cell_info_parsed_df[,"cell_type"]=="S") & (cell_info_parsed_df[,"hdr_info"]%in%c("WT", "HDR")) ),]

cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"]=as.character(cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"])
cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"]=gsub("WT", "Human", cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"])
cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"]=gsub("HDR", "Chimp-Edited", cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"])
cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"]=factor(cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"], c("Chimp-Edited", "Human"))

plotName=paste("UMAP_LOXL2_exp_WT_HDR_S_only", sep="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_info_parsed_df_s_wt_hdr_only, aes_string(x="UMAP_1", y="UMAP_2", colour="exp", stroke="hdr_info"))+geom_point(alpha=0.5, size=7)+scale_colour_gradient(name="LOXL2 Expression", low = "grey85", high = "red4", aesthetics = "colour")+scale_discrete_manual(aesthetics = "stroke", name="", values = c("Human" = 1, "Chimp-Edited" = 3))+xlab("UMAP 1")+ylab("UMAP 2")+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
ggsave(filename=plotOutFileName, plot=p, width=20, height=20)


##### make plots of LOXL2 expression  #####
plotName=paste("violin_LOXL2_exp_WT_HDR_S_only", sep="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_info_parsed_df_s_wt_hdr_only, aes_string(x="hdr_info", y="exp", stroke="hdr_info"))+geom_violin()+geom_point(alpha=0.5, size=7, position="jitter")+geom_signif(comparisons = list(c("Chimp-Edited", "Human")), test="wilcox.test", textsize=8, map_signif_level=FALSE, tip_length=0, y_position=3)+scale_discrete_manual(aesthetics = "stroke", name="", values = c("Human" = 1, "Chimp-Edited" = 5))+xlab("Type")+ylab("LOXL2 Expression")+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
ggsave(filename=plotOutFileName, plot=p, width=20, height=20)


plotName=paste("UMAP_LOXL2_exp_WT_HDR_S_only_alt", sep="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_info_parsed_df_s_wt_hdr_only, aes_string(x="UMAP_1", y="UMAP_2", colour="hdr_info"))+geom_point(alpha=0.5, size=7)+scale_colour_discrete(name="HDR Info")+xlab("UMAP 1")+ylab("UMAP 2")+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
ggsave(filename=plotOutFileName, plot=p, width=20, height=20)

######### address review concern on subsampling #########

color_mapping=c("#E5AC23", "#1B335F")
names(color_mapping)=c("Human", "Chimp-Edited")

plotName=paste("density_LOXL2_exp_WT_HDR_S_only_alt", sep="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_info_parsed_df_s_wt_hdr_only, aes_string(x="UMAP_2", y="UMAP_1"))+stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = hdr_info), bins=5)+xlab("UMAP 2")+ylab("UMAP 1")+scale_fill_manual(values=color_mapping)+guides(colour= guide_legend(title = ""))+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
ggsave(filename=plotOutFileName, plot=p, width=20, height=20)

table(as.character(cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"]))

#points with density, fig. S10C
plotName=paste("UMAP_with_density_LOXL2_exp_WT_HDR_S_only_alt", sep="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
p=ggplot(cell_info_parsed_df_s_wt_hdr_only, aes_string(x="UMAP_2", y="UMAP_1", colour="hdr_info"))+geom_point(alpha=0.5, size=7)+scale_colour_manual(name="", values=color_mapping)+stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = hdr_info), bins=5)+xlab("UMAP 2")+ylab("UMAP 1")+scale_fill_manual(name="", values=color_mapping)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
ggsave(filename=plotOutFileName, plot=p, width=20, height=20)

#plot using actual subsamples
color_mapping=c("#E5AC23", "#1B335F")
names(color_mapping)=c("Human", "Chimp-Edited")

num_subsamples=3

for(i in 1:num_subsamples){
	chimp_ind=which(as.character(cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"])=="Chimp-Edited")
	human_ind=which(as.character(cell_info_parsed_df_s_wt_hdr_only[,"hdr_info"])=="Human")
	
	#subsample human ind
	set.seed(300+i)
	rand_sample_human_ind=sample(human_ind, length(chimp_ind))
	
	all_ind=c(rand_sample_human_ind, chimp_ind)
	
	plotName=paste("Random_Sample", i, "UMAP_LOXL2_exp_WT_HDR_S_only_alt", sep="_")
	plotTitle=plotName
	plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
	p=ggplot(cell_info_parsed_df_s_wt_hdr_only[all_ind,], aes_string(x="UMAP_2", y="UMAP_1", colour="hdr_info"))+geom_point(alpha=0.5, size=7)+xlab("UMAP 2")+ylab("UMAP 1")+scale_colour_manual(name="HDR Info", values=color_mapping)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=40), axis.text.x=element_text(size=40), axis.title.y=element_text(size=40), axis.text.y=element_text(size=40), legend.title=element_text(size=40), legend.text=element_text(size=40), legend.key.size = unit(1.5, "cm"), legend.key = element_rect(colour = NA, fill = NA), axis.line.x = element_line(color="black", size=1), axis.line.y = element_line(color="black", size=1))
	ggsave(filename=plotOutFileName, plot=p, width=20, height=20)
}

#******************** differential expression via DESeq2, volcano plot/bar plots ********************#

wt_ind=which((cell_info@meta.data$hdr_info=="WT") & (cell_info@meta.data$cell_type=="S"))
wt_ids=rownames(cell_info@meta.data)[wt_ind]

hdr_ind=which((cell_info@meta.data$hdr_info=="HDR") & (cell_info@meta.data$cell_type=="S"))
hdr_ids=rownames(cell_info@meta.data)[hdr_ind]

#run DESEQ2 after above # 

test_use="DESeq2"
cluster1.markers <- FindMarkers(cell_info, ident.1 = wt_ids, ident.2 = hdr_ids, logfc.threshold=0.1, test.use=test_use)
head(cluster1.markers, 100)

cluster1.markers.out=cluster1.markers
clust_gene_names=rownames(cluster1.markers.out)
cluster1.markers.out$gene=clust_gene_names
rownames(cluster1.markers.out)=c()
cluster1.markers.out=cluster1.markers.out[,c("gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")]

file_out_name=paste(out_file_full_path, "_", "diff_markers_wt_hdr", "_", test_use, ".txt", sep="")

write.table(cluster1.markers.out, file_out_name, quote=F, row.names=F)

#Fig. 4F volcano plot
padj_cutoff=0.1

file_name="/home_path/10X/Analysis_Output/del_1a_2_13_21_diff_markers_wt_hdr_DESeq2.txt"
deseq_results=read.table(file_name, header=T)
deseq_results[,"avg_logFC_2"]=log(exp(deseq_results[,"avg_logFC"]), base=2)
length(which(deseq_results[,"p_val_adj"]<padj_cutoff))
select_labels=as.character(deseq_results[which(deseq_results[,"p_val_adj"]<padj_cutoff),"gene"])

plotName=paste(c("hdr_vs_unedited_deseq2_volcano", "padj_cutoff", padj_cutoff), collapse="_")
plotTitle=plotName
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
pdf(plotOutFileName, onefile=FALSE)
EnhancedVolcano(deseq_results, lab=as.character(deseq_results[,"gene"]), selectLab=select_labels, x="avg_logFC_2", y="p_val_adj", title = 'HDR (Chimp) vs. Unedited (Human)', subtitle = "", caption = paste0('Total = ', nrow(deseq_results), ' Genes'), pCutoff=0.1, FCcutoff=5)
dev.off()

######## go enrichment, bar plots ##########

gene_diff_vec_keep=rep(0, nrow(deseq_results))
names(gene_diff_vec_keep)=as.character(deseq_results[,"gene"])

diff_genes_ind=which(deseq_results[,"p_val_adj"]<0.1)
gene_diff_vec_keep[diff_genes_ind]=1

pwf=nullp(gene_diff_vec_keep,"hg19","geneSymbol")

GO.nobias=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
GO.nobias[,"FDR"]=p.adjust(GO.nobias$over_represented_pvalue, method="BH")

enriched.GO=GO.nobias[which(GO.nobias[,"FDR"]<.05),]
enriched.GO[,"log_FDR"]=-log(enriched.GO[,"FDR"], base=10)
enriched.GO=enriched.GO[order(enriched.GO[,"log_FDR"]),]
enriched.GO[,"term"]=factor(as.character(enriched.GO[,"term"]), as.character(enriched.GO[,"term"]))

development_ind=which(grepl("neur|development", enriched.GO[,"term"]))
development.enriched.GO=enriched.GO[development_ind,]

motility_ind=which(grepl("adhesion|migration|motility", enriched.GO[,"term"]))
motility.enriched.GO=enriched.GO[motility_ind,]

p=ggplot(development.enriched.GO, aes(x = factor(term), y=log_FDR)) +geom_bar(stat="identity")+xlab("")+ylab("log FDR")+ theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))+coord_flip()
plotName=paste("development_enriched_GO_FDR", sep="_")
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
ggsave(filename=plotOutFileName, plot=p, width=10, height=10 )

p=ggplot(motility.enriched.GO, aes(x = factor(term), y=log_FDR)) +geom_bar(stat="identity")+xlab("")+ylab("log FDR")+ theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))+coord_flip()
plotName=paste("motility_enriched_GO_FDR", sep="_")
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
ggsave(filename=plotOutFileName, plot=p, width=10, height=10 )

combined.enriched.GO=rbind(development.enriched.GO, motility.enriched.GO)

combined.enriched.GO[,"term"]=factor(as.character(combined.enriched.GO[,"term"]), as.character(combined.enriched.GO[,"term"]))
p=ggplot(combined.enriched.GO, aes(x = factor(term), y=log_FDR)) +geom_bar(stat="identity")+xlab("")+ylab("log FDR")+ theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))+coord_flip()
plotName=paste("combined_enriched_GO_FDR", sep="_")
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
ggsave(filename=plotOutFileName, plot=p, width=20, height=10 )

#Fig. 4D
combined.enriched.GO[,"term"]=factor(as.character(combined.enriched.GO[,"term"]), as.character(combined.enriched.GO[,"term"]))
p=ggplot(combined.enriched.GO, aes(x = factor(term), y=numDEInCat, fill=FDR)) +geom_bar(stat="identity")+xlab("")+ylab("Number of Genes")+ theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+coord_flip()
plotName=paste("combined_enriched_GO_FDR_num_genes", sep="_")
plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
ggsave(filename=plotOutFileName, plot=p, width=20, height=10 )


