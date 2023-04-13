#this script takes a bed intersected file of labels, and draws a random label

import sys
import random

arg = sys.argv

intersectedFileName=arg[1]
outFileName=arg[2]

outFile=open(outFileName, "w")

#read in coordinates
prevDelCoordId=""
with open(intersectedFileName, "r") as intersectedFile:
	for ind, l1 in enumerate(intersectedFile):
					 
		l1=l1.split()
		chrName=l1[0]
		startPos=int(l1[1])
		endPos=int(l1[2])
		delCoordId=l1[3]
		
		labelChr=l1[4]
		labelStartPos=int(l1[5])
		endStartPos=int(l1[6])
		
		labelCategory=l1[7]
		
		if labelChr==".":
			outFile.write(delCoordId+"\t"+"NA"+"\n")
			continue
		#first coordinate with overlapped results
		if prevDelCoordId=="":
			prevDelCoordId=delCoordId
			coordLabels=[labelCategory]
		else:	
			#switched labels	
			if prevDelCoordId!=delCoordId:
				#write out previous label
				labelKeep=coordLabels[random.randint(0,len(coordLabels)-1)]
				outFile.write(prevDelCoordId+"\t"+labelKeep+"\n")
				#reset
				prevDelCoordId=delCoordId
				coordLabels=[labelCategory]
			#prevDelCoordId==delCoordId
			else:
				coordLabels.append(labelCategory)

#write out last position
labelKeep=coordLabels[random.randint(0,len(coordLabels)-1)]
outFile.write(prevDelCoordId+"\t"+labelKeep+"\n")

outFile.close()

 
