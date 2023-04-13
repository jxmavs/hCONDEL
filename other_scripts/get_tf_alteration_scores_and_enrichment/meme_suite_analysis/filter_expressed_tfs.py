import sys
import numpy as np

arg=sys.argv

expressedTFFileName=arg[1]
TFDBFileName=arg[2]
outFileName=arg[3]

expressedTFDict=dict()
with open(expressedTFFileName, "r") as expressedTFFile:
	for ind, l1 in enumerate(expressedTFFile):
		l1=l1.rstrip().split("\t")
		gene_ensembl_id=l1[0]
		gene_name=l1[1]
		expressedTFDict[gene_name]=gene_ensembl_id

#tfName is tf name from motif database
#tfNameMatchEnsembl is tf name that should match the gene name from ensembl 

tfArr=[]
with open(TFDBFileName, "r") as TFDBFile:
	for ind, l1 in enumerate(TFDBFile):
		l1=l1.rstrip().split("\t")
		tfID=l1[0]
		tfName=l1[1]
		tfNameMatchEnsembl=l1[2]
		tfArr.append([tfID, tfName, tfNameMatchEnsembl])

outFile=open(outFileName, "w")
for ind in range(len(tfArr)):
	tfID=tfArr[ind][0]
	tfName=tfArr[ind][1]
	tfNameMatchEnsembl=tfArr[ind][2]
	expressed=1
	#if TF complex, then check that each TF is expressed
	if("::" in tfNameMatchEnsembl): 
			
		tfComplexArr=tfNameMatchEnsembl.split("::")
		for ind in range(len(tfComplexArr)):
			if(tfComplexArr[ind] not in expressedTFDict):
				expressed=0
				break
	else:
		if(tfNameMatchEnsembl not in expressedTFDict):
			expressed=0
	
	if(expressed==1):
		outFile.write(tfName+"\n")

outFile.close()
