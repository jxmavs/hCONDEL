#remove human annotated deleted sequence from fasta

import sys

arg=sys.argv

fastaFileName=arg[1]
outFileName=arg[2]
idInfoAdd=arg[3]

outFile=open(outFileName, "w")

with open(fastaFileName, "r") as fastaFile:
	for ind, l1 in enumerate(fastaFile):
		l1=l1.split()[0]
			
		if ind%2==0:
			coordIdInfo=l1
			coordIdInfo=coordIdInfo.split("#")
			origId=coordIdInfo[0]
			origId=origId+"_"+idInfoAdd
			coordIdInfo[0]=origId
			coordIdInfoNew="#".join(coordIdInfo)	
		else:
			sequence=l1	
			
			outFile.write(coordIdInfoNew+"\n"+sequence+"\n")

outFile.close()

