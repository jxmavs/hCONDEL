#this file is cleaned from the original version in MPRADel_NovaSeqAnalysis_v3.r
#all_cell_types_combined_mpradel_deseq.txt includes the positive controls from Tewhey et al. 2016, which includes ref/alt variants
#here ref is chimp, alt is human

library(DESeq2) #Differential Expression Analysis

#*********************************** always run the beginning part here in every analysis for other analyses below ************************************#

#columns are:  treat type, ref type, comparison group, row subset group, name of comparison group 

inputConditionsFileName="/home_path/all_cell_types_combined_mpradel_deseq_cond_file.txt"
inputCountFileName="/home_path/all_cell_types_combined_mpradel_deseq.txt"

compareInfoDFFileName="/home_path/NovaSeq_Analysis/novaseq_cell_type_comparisons_v4.txt"

setwd("/home_path/NovaSeq_Analysis")
out_file_path="/home_path/NovaSeq_Analysis/Analysis_Output"
out_prefix="novaseq_12_15_18"

out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")

##########################################################################

#MPRA DESeq Analysis

countData=read.table(inputCountFileName, header=TRUE, row.names=1, check.names=FALSE)
colData=read.table(inputConditionsFileName, header=FALSE, row.names=1)
compareInfoDF=read.table(compareInfoDFFileName, header=FALSE)
colnames(colData)="CellType"


#for DESEQ for finding ref/alt oligos
ref_tag="CC$|ref$"
alt_tag="HH$|alt$"

################# initial processing of files #############

rowNames=rownames(countData)
sortedCountDataInd=order(rowNames)

countDataParsed=countData[sortedCountDataInd,]


###################################################

#have a seperate countData for each comparison 
#count data should be sorted with ref first/alt second

numComparisonSets=max(compareInfoDF[,3])
countDataRowNames=rownames(countDataParsed)

countDataAllSets=list()
colDataAllSets=list()

for(ind in 1:numComparisonSets){
	
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	
	treat_and_nontreat_cols=unique(c(as.character(compareInfoDFSubset[,1]),as.character(compareInfoDFSubset[,2])))
	
	#get col data names
	col_subset_ind=rownames(colData)[which(colData[,1] %in% treat_and_nontreat_cols)]
	
	countDataSubset=countDataParsed[, col_subset_ind]
	colDataSubset=colData[which(colData[,1]%in%treat_and_nontreat_cols),,drop=FALSE]
	
	countDataAllSets[[ind]]=countDataSubset
	colDataAllSets[[ind]]=colDataSubset
}

#*****************************************  using DESEQ2 to calculate FC/Skew ********************************************************#
allDESEQResults=list()

#boolean to make plots or not
plotDESEQ=TRUE

for(ind in 1:numComparisonSets){
	countDataSubset=countDataAllSets[[ind]]
	compareInfoDFSubset=compareInfoDF[which(compareInfoDF[,3]==ind),]
	colDataSubset=colDataAllSets[[ind]]
	
	ref_ind=which(grepl(ref_tag, rownames(countDataSubset)))
	alt_ind=which(grepl(alt_tag, rownames(countDataSubset)))
	
	Ref_Count_Data=countDataSubset[ref_ind,]
	Alt_Count_Data=countDataSubset[alt_ind,]
	
	colnames(Ref_Count_Data)=paste(colnames(Ref_Count_Data), "Ref", sep="_")
	colnames(Alt_Count_Data)=paste(colnames(Alt_Count_Data), "Alt", sep="_")

	refRowNames=rownames(Ref_Count_Data)
	refID=strsplit(as.character(refRowNames),  "\\|")
	refSeqNames=sapply(refID, function(x){if( ("EMVAR" %in% x) | ("LCL_HEPG2_POS" %in% x) ){ paste(x[4:(length(x)-1)], collapse="|") } else{ paste(x[7:(length(x)-1)], collapse="|") } } )

	#### sanity check that ref/alt row names are same ####
	altRowNames=rownames(Alt_Count_Data)
	altID=strsplit(as.character(refRowNames),  "\\|")
	altSeqNames=sapply(altID, function(x){if( ("EMVAR" %in% x) | ("LCL_HEPG2_POS" %in% x) ){ paste(x[4:(length(x)-1)], collapse="|") } else{ paste(x[7:(length(x)-1)], collapse="|") } } )

	if(!all(refSeqNames==altSeqNames)){
		print("ref and alt names don't match")
		break
	}
	
	#########################################################
	
	rownames(Ref_Count_Data)=refSeqNames
	rownames(Alt_Count_Data)=refSeqNames

	#combine columns
	revisedCountData=data.frame(Ref_Count_Data, Alt_Count_Data)

	#revise colData, add Ref/Alt info 
	revisedColData=rbind(colDataSubset, colDataSubset)
	rownames(revisedColData)[1:(nrow(revisedColData)/2)]=colnames(Ref_Count_Data)
	rownames(revisedColData)[(nrow(revisedColData)/2+1):nrow(revisedColData)]=colnames(Alt_Count_Data)

	revisedColData$Ref_Alt=c( rep("Ref", nrow(revisedColData)/2 ), rep("Alt", nrow(revisedColData)/2) )

	revisedColData[,"Ref_Alt"]=factor(revisedColData[,"Ref_Alt"])
	revisedColData[,"CellType"]=factor(revisedColData[,"CellType"])

	#loop and run DESEQ across all cell type comparisons 
	deseqDataTableRevisedAll=data.frame()

	for(i in 1:nrow(compareInfoDFSubset)){

		treat_level=as.character(compareInfoDFSubset[i,1])
		base_level=as.character(compareInfoDFSubset[i,2])
		
		revisedColData[,"Ref_Alt"]=relevel(revisedColData[,"Ref_Alt"], "Ref")
		revisedColData[,"CellType"]=relevel(revisedColData[,"CellType"], base_level)

		deseqData=DESeqDataSetFromMatrix(countData = revisedCountData, colData = revisedColData, design = ~ Ref_Alt+CellType+Ref_Alt:CellType)
		#model.matrix(design(deseqData), colData(deseqData))

		deseqDataResults=deseqData
		sizeFactors(deseqDataResults)=estimateSizeFactorsForMatrix(revisedCountData)
				
		deseqDataResults<-nbinomWaldTest(deseqDataResults, betaPrior=FALSE)

		resultsNames(deseqDataResults)
		
		deseqOutputSkew=results(deseqDataResults, name=paste(c("Ref_AltAlt.CellType", treat_level), collapse=""))
		
		deseqOutputRef=results(deseqDataResults, contrast=c("CellType", treat_level, base_level))
		
		deseqOutputAlt=results(deseqDataResults, contrast=list(c(paste(c("CellType", treat_level, "vs", base_level), collapse="_"), paste(c("Ref_AltAlt.CellType", treat_level), collapse="") )) )

		#Skew

		deseqOutputDataFrame=data.frame(deseqOutputSkew)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedSkew=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedSkew)[2:ncol(deseqDataTableRevisedSkew)]=paste(colnames(deseqDataTableRevisedSkew)[2:ncol(deseqDataTableRevisedSkew)], "Skew", treat_level, sep="_")

		#Ref 

		deseqOutputDataFrame=data.frame(deseqOutputRef)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedRef=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedRef)[2:ncol(deseqDataTableRevisedRef)]=paste(colnames(deseqDataTableRevisedRef)[2:ncol(deseqDataTableRevisedRef)], "Ref", treat_level, sep="_")

		#Alt

		deseqOutputDataFrame=data.frame(deseqOutputAlt)
		OriginalName=rownames(deseqOutputDataFrame)
		rownames(deseqOutputDataFrame)=c()

		deseqDataTableRevisedAlt=data.frame(SeqName=OriginalName, deseqOutputDataFrame)
		colnames(deseqDataTableRevisedAlt)[2:ncol(deseqDataTableRevisedAlt)]=paste(colnames(deseqDataTableRevisedAlt)[2:ncol(deseqDataTableRevisedAlt)], "Alt", treat_level, sep="_")

		deseqDataTableRevised_1=merge(deseqDataTableRevisedSkew, deseqDataTableRevisedRef, by="SeqName")
		deseqDataTableRevised=merge(deseqDataTableRevised_1, deseqDataTableRevisedAlt, by="SeqName")

		if(nrow(deseqDataTableRevisedAll)==0){
			deseqDataTableRevisedAll=deseqDataTableRevised
		}
		else{
			deseqDataTableRevisedAll=merge(deseqDataTableRevisedAll, deseqDataTableRevised, by="SeqName")
		}
	
		####################### some prelim DESEQ plots ##########################
		
		if(plotDESEQ){
	
			plotName=paste(c("treat", treat_level, "base", base_level,"hCONDEL_Skew"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputSkew, alpha=0.01, ylab="Skew", xlab="Mean of Normalized Counts")
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"hCONDEL_Ref_Activity"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputRef, alpha=0.01, ylab="Ref (Chimp) Activity", xlab="Mean of Normalized Counts" )
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"hCONDEL_Alt_Activity"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotMA(deseqOutputAlt, alpha=0.01, ylab="Alt (Human) Activity", xlab="Mean of Normalized Counts" )
			dev.off()

			plotName=paste(c("treat", treat_level, "base", base_level,"disp_fit"), collapse="_")
			plotTitle=plotName
			plotOutFileName=paste(c(out_file_full_path,paste(c("_",plotName),collapse=""), ".pdf"), collapse="")
			pdf(plotOutFileName)
			par(cex.axis=1, cex.lab=1.5, cex.main=1, cex.sub=1)
			plotDispEsts(deseqDataResults)
			dev.off()
			
		}
	}
	
	allDESEQResults[[ind]]=deseqDataTableRevisedAll

}

### combine table into one file ###

deseq_df_all=allDESEQResults[[1]]
colnames_deseq=colnames(deseq_df_all)
#keep padj and log2FoldChange
keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*|baseMean*", colnames_deseq))
deseq_df_all=deseq_df_all[,c(1, keep_col_ind)]

#collapse to single data frame
for(ind in 2:length(allDESEQResults)){
	
	deseq_df_iter=allDESEQResults[[ind]]
	colnames_deseq=colnames(deseq_df_iter)
	#keep padj and log2FoldChange
	keep_col_ind=which(grepl("padj*|log2FoldChange*|lfcSE*|pvalue*|baseMean*", colnames_deseq))
	
	deseq_df_all=merge(deseq_df_all, deseq_df_iter[, c(1, keep_col_ind)], by="SeqName")
}

#write out deseq table
out_file_full_path=paste(c(out_file_path, out_prefix), collapse="/")
textOutFileName=paste(c(out_file_full_path, "_deseq_results_all_table.txt"), collapse="")
write.table(deseq_df_all, quote=F, row.names=F, textOutFileName)


