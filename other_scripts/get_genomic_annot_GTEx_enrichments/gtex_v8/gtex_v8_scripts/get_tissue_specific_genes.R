library(data.table)
library(limma)
library(edgeR)
args = commandArgs(trailingOnly=TRUE)

sample_gene_exp_table_file_name=args[1]
sample_attributes_table_file_name=args[2]
out_file_path_with_file_header=args[3]
tissue_of_interest=args[4]

##### pre-install packages if necessary ####
#change path in ~/.Renviron to /cluster_path/bin/libs/Rlibs/R-4.0
#wget -P /cluster_path/bin/libs/Rlibs/R-4.0 https://bioconductor.org/packages/release/bioc/src/contrib/limma_3.52.1.tar.gz
#R CMD INSTALL -l /cluster_path/bin/libs/Rlibs/R-4.0 /cluster_path/bin/libs/Rlibs/R-4.0/limma_3.52.1.tar.gz

#wget -P /cluster_path/bin/libs/Rlibs/R-4.0 https://bioconductor.org/packages/release/bioc/src/contrib/edgeR_3.38.1.tar.gz
#R CMD INSTALL -l /cluster_path/bin/libs/Rlibs/R-4.0 /cluster_path/bin/libs/Rlibs/R-4.0/edgeR_3.38.1.tar.gz

##### for testing #####
#sample_gene_exp_table_file_name="/cluster_path/ape_project/deletions_project/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
#head -1000 /cluster_path_temp/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct> /cluster_path_temp/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.1000
#sample_gene_exp_table_file_name="/cluster_path_temp/gtex_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.1000"
#sample_attributes_table_file_name="/cluster_path/ape_project/deletions_project/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.subset.count_table_intersect.txt"
#out_file_path_with_file_header="/cluster_path/ape_project/deletions_project/gtex_v8/output/novaseq_12_15_18_denylist_filtered_broad_tissues"
#tissue_of_interest="Brain"

#use for gzipped file, too slow
#sample_gene_exp_table_file_gz=gzfile(sample_gene_exp_table_file_name,'rt')  
#gene_exp_table=read.table(sample_gene_exp_table_file_gz, header=T, skip=2)

gene_exp_table=as.data.frame(fread(sample_gene_exp_table_file_name, header=T, skip=2, sep="\t"))

sample_attributes_table=read.table(sample_attributes_table_file_name, header=T, sep="\t", quote="")
colnames(sample_attributes_table)=c("sample_id", "broad_tissue", "specific_tissue")

#get tissue map
tissue_type_map=unique(sample_attributes_table[,c("broad_tissue", "specific_tissue")])

#get tissue types
all_broad_tissues=unique(as.character(sample_attributes_table[,"broad_tissue"]))
all_specific_tissues=unique(as.character(sample_attributes_table[,"specific_tissue"]))

#check to see if tissue is specific of broad
if(tissue_of_interest %in% all_broad_tissues ){
	analysis_type="broad"
} else if (tissue_of_interest %in% all_specific_tissues ){
	analysis_type="specific"
} else{
	print("Invalid Tissue")
	print(tissue_of_interest)
	quit()
}

gene_exp_table_gene_ids_table=gene_exp_table[,c(1,2)]
colnames(gene_exp_table_gene_ids_table)=c("ensembl_id", "gene_name")

if(analysis_type=="broad"){
	
	#get samples corresponding to tissue of interest
	treatment_ind=which(sample_attributes_table[,"broad_tissue"]==tissue_of_interest)
	treatment_sample_ids=as.character(sample_attributes_table[treatment_ind, "sample_id"])

	#get sample ids from other tissues
	control_ind=which(sample_attributes_table[,"broad_tissue"]!=tissue_of_interest)
	control_sample_ids=as.character(sample_attributes_table[control_ind, "sample_id"])
	
} else{
	
	#specific tissue
	#get samples corresponding to tissue of interest
	treatment_ind=which(sample_attributes_table[,"specific_tissue"]==tissue_of_interest)
	treatment_sample_ids=as.character(sample_attributes_table[treatment_ind, "sample_id"])
	
	broad_tissue=as.character(tissue_type_map[which(tissue_type_map[,"specific_tissue"]==tissue_of_interest), "broad_tissue"])
	
	if(tissue_of_interest=="Brain - Cortex"){
		#get sample ids from other specific tissues in the same class as the tissue of interest
		all_other_tissue_ind=which((sample_attributes_table[,"broad_tissue"]==broad_tissue) & (sample_attributes_table[,"specific_tissue"]!=tissue_of_interest) & (sample_attributes_table[,"specific_tissue"]!="Brain - Frontal Cortex (BA9)"))
		control_sample_ids=as.character(sample_attributes_table[all_other_tissue_ind, "sample_id"])
			
	} else if(tissue_of_interest=="Brain - Frontal Cortex (BA9)"){
		#get sample ids from other specific tissues in the same class as the tissue of interest
		all_other_tissue_ind=which((sample_attributes_table[,"broad_tissue"]==broad_tissue) & (sample_attributes_table[,"specific_tissue"]!=tissue_of_interest) & (sample_attributes_table[,"specific_tissue"]!="Brain - Cortex"))
		control_sample_ids=as.character(sample_attributes_table[all_other_tissue_ind, "sample_id"])
		
	} else if(tissue_of_interest=="Brain - Cerebellum"){
        #get sample ids from other specific tissues in the same class as the tissue of interest
        all_other_tissue_ind=which((sample_attributes_table[,"broad_tissue"]==broad_tissue) & (sample_attributes_table[,"specific_tissue"]!=tissue_of_interest) & (sample_attributes_table[,"specific_tissue"]!="Brain - Cerebellar Hemisphere"))
        control_sample_ids=as.character(sample_attributes_table[all_other_tissue_ind, "sample_id"])

    } else if(tissue_of_interest=="Brain - Cerebellar Hemisphere"){
        #get sample ids from other specific tissues in the same class as the tissue of interest
        all_other_tissue_ind=which((sample_attributes_table[,"broad_tissue"]==broad_tissue) & (sample_attributes_table[,"specific_tissue"]!=tissue_of_interest) & (sample_attributes_table[,"specific_tissue"]!="Brain - Cerebellum"))
        control_sample_ids=as.character(sample_attributes_table[all_other_tissue_ind, "sample_id"])

    } else{
		#get sample ids from other specific tissues in the same class as the tissue of interest
		all_other_tissue_ind=which((sample_attributes_table[,"broad_tissue"]==broad_tissue) & (sample_attributes_table[,"specific_tissue"]!=tissue_of_interest))
		control_sample_ids=as.character(sample_attributes_table[all_other_tissue_ind, "sample_id"])
	}
}

######## preprocess #############
#adapted from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4
gene_exp_table_subset=gene_exp_table[,c(treatment_sample_ids, control_sample_ids)]
rownames(gene_exp_table_subset)=as.character(gene_exp_table[,1])
conditions<-factor(c(rep("treatment", length(treatment_sample_ids)), rep("control", length(control_sample_ids))))

## get cpms ##
y <- DGEList(counts=gene_exp_table_subset, group=conditions)
##Remove rows consistently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

## get p-values, FDR ##
pvalues <- sapply(1:nrow(count_norm),function(i){
     data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
     p=wilcox.test(gene~conditions, data)$p.value
     return(p)
   })
fdr=p.adjust(pvalues,method = "fdr")

## get l2fcs ##
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

## create out_table ##

#print out number of kept genes
print(length(which(keep)))
#sanity check that gene names are the same as row names
print("gene name check")
print(all(rownames(count_norm)==as.character(gene_exp_table_gene_ids_table[keep,1])))

out_table_temp_1=data.frame(gene_exp_table_gene_ids_table[keep,], log2FC=foldChanges, p_val=pvalues, FDR=fdr)
#add on genes with low counts that were filtered out if it exists, put NAs in other columns
if(length(which(!keep))!=0){
	out_table_temp_2=data.frame(gene_exp_table_gene_ids_table[!keep,], log2FC=NA, p_val=NA, FDR=NA)
	out_table=rbind(out_table_temp_1, out_table_temp_2)
} else{
	out_table=out_table_temp_1
}
out_table=out_table[, c("ensembl_id", "gene_name", "log2FC", "p_val", "FDR")]
out_table=out_table[order(as.character(out_table[,"ensembl_id"])),]

#write out gene names that passed FDR correction
textOutFileName=paste(paste(c(out_file_path_with_file_header, gsub(" ", "_", tissue_of_interest)), collapse="_"), ".txt", sep="")
write.table(out_table, textOutFileName, quote=F, row.names=F, sep="\t")

print(warnings())




