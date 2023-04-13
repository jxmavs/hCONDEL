import sys

#this script gets the average signal score/ sees if any deleitons intersect a region of interest
arg = sys.argv

intersectedFileName=arg[1]
cutOffVal=float(arg[2])
outFileName=arg[3]

intersectDict=dict()

with open(intersectedFileName, "r") as intersectedFile:
	for ind, l1 in enumerate(intersectedFile):
		l1=l1.split()
		coordID='\t'.join(l1[0:3])
		#print(l1)
		chrom=l1[4]
		startPos=int(l1[5])
		endPos=int(l1[6])
		score=l1[len(l1)-1]
		#print(coordID) 
		if(coordID not in intersectDict):
			#first array is the coordinates that intsersect, second array is the score
			intersectDict[coordID]=[[],[]]
			if(startPos!=-1):		
				intersectDict[coordID][0].append([chrom, startPos, endPos, float(score)])
		#the ID is already in the dictionary, so simply append on more info
		else:
			intersectDict[coordID][0].append([chrom, startPos, endPos, float(score)])

#now iterate through intersectDict and calculate average score as well as checking to see if any overlap
#allKeys=intersectDict.keys()
for key in intersectDict:
	coordArrs=intersectDict[key][0]
	#skip calculations of any coordinates which don't overlap with anything
	#output -1, -1, -1
	if(coordArrs==[]):
		intersectDict[key][1]=[-1,-1,-1, -1]
		continue
	totalLen=0
	totalScore=0
	containsSigScore=0
	#totalNumCoords is the number of intervals that intersect 
	totalNumCoords=len(coordArrs)
	for ind in range(totalNumCoords):
		coord=coordArrs[ind]
		startPos=coord[1]
		endPos=coord[2]
		score=coord[len(coord)-1]
		length=endPos-startPos
		
		totalScore=length*score+totalScore
		totalLen=length+totalLen
		
		if(score>=cutOffVal):
			containsSigScore=1
	#avgScore=totalScore/float(totalLen)
	intersectDict[key][1]=[totalScore, totalLen, containsSigScore, totalNumCoords]

with open(outFileName, "w") as outFile:      
	for key in intersectDict:
		outFile.write(key+"\t")
		
		scoreInfo=intersectDict[key][1]
		for ind in range(len(scoreInfo)):
			outFile.write(str(scoreInfo[ind])+"\t")
		outFile.write("\n")
