#this script combines output from closest bed so taht same distance entries are merged
#this script assumes that the second column contains the distance, and all subsequent columns after are meta columns

import sys

arg=sys.argv

closestBedFileName=arg[1]
outFileName=arg[2]

outFileInfo=dict() 
with open(closestBedFileName, "r") as closestBedFile:
	for ind, l1 in enumerate(closestBedFile):
		l1=l1.split()
		idInterest=l1[0]
		otherCols=l1[1:len(l1)]
			
		#if multiple entries with same distance, then merge
		if (idInterest in outFileInfo):
			#combine fields separated by "|"
			for ind in range(len(outFileInfo[idInterest])):
				newField=outFileInfo[idInterest][ind]+"|"+otherCols[ind]
				outFileInfo[idInterest][ind]=newField
		else:
			outFileInfo[idInterest]=otherCols

with open(outFileName, "w") as outFile:
	for idInterest in outFileInfo:
		outFile.write(idInterest+"\t")
		for ind in range(len(outFileInfo[idInterest])):
			outFile.write(str(outFileInfo[idInterest][ind])+"\t")
		outFile.write("\n")


