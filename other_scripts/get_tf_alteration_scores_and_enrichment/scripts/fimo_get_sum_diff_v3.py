#this script parses FIMO output and obtains the max difference between comparison files
#DON'T USE ${scriptdir}/fimo_get_sum_diff.py as there is a typo with calculating centerStartCoord and centerEndCoord
#and the deletion overlap only had the situation of complete tf overlapping del
import sys
import numpy as np

arg=sys.argv

fimoInputFileName=arg[1]
pairwiseCombFileName=arg[2]
outFileName=arg[3]
lookAtCenteredOnly=int(arg[4])
cutOffVal=float(arg[5])

allSpecies=[]
pairwiseArr=[]
with open(pairwiseCombFileName, "r") as pairwiseCombFile:
    for ind, l1 in enumerate(pairwiseCombFile):
        l1=l1.split()
        pairwiseArr.append(l1)
        #hacky way to obtain all unique species going to be used in pairwise comparison
        allSpecies.append(l1[0])
        allSpecies.append(l1[1])

allSpecies=set(allSpecies)

allTFDict=dict()
with open(fimoInputFileName, "r") as fimoInputFile:
	for ind, l1 in enumerate(fimoInputFile):
		l1=l1.split()
		tfID=l1[0]
		coordIdInfo=l1[1]
		startPos=int(l1[2])
		endPos=int(l1[3])
		strand=l1[4]
		score=float(l1[5])
		pval=float(l1[6])
		
		
		#I can check to see if the motif lies exactly over the deletion position if necessary
		         
		coordIdInfo=coordIdInfo.split("#")
		speciesInfo=coordIdInfo[1]
		speciesInfo=speciesInfo.split("|")
		speciesID=speciesInfo[0]
		
		#remember that currently the position is 1 bp after the deletion, and speciesInfo[5] is equivalent to speciesInfo[6] for humans
		if speciesID=="hg38":
			#bp before deletion
			centerStartCoord=int(speciesInfo[5])-int(speciesInfo[2])
			#bp after deletion
			centerEndCoord=int(speciesInfo[6])-int(speciesInfo[2])+1
		else:
			centerStartCoord=int(speciesInfo[5])-int(speciesInfo[2])+1
			centerEndCoord=int(speciesInfo[6])-int(speciesInfo[2])
		
		#if we have the lookAtCenteredOnly==1 option set, then only look at the max centered on the deletion site      
		if(lookAtCenteredOnly==1):
			if speciesID=="hg38":
				if( (endPos<=centerStartCoord) or (centerEndCoord<=startPos) ):
					continue
			else:
				if( (endPos<centerStartCoord) or (centerEndCoord<startPos) ):
					continue
		
		#if( (startPos<=centerStartCoord) and (centerEndCoord<=endPos)):
		#    overlapDel=1
		#else:
		#    overlapDel=0		
		
		coordId=coordIdInfo[0]
		
		if(coordId not in allTFDict):
			allTFDict[coordId]=dict()
		if(tfID not in allTFDict[coordId]):
			allTFDict[coordId][tfID]=dict()
			allTFDict[coordId][tfID][speciesID]=[]
			if(pval<cutOffVal):
				allTFDict[coordId][tfID][speciesID].append(score)
		elif(speciesID not in allTFDict[coordId][tfID]):
			allTFDict[coordId][tfID][speciesID]=[]
			if(pval<cutOffVal):
				allTFDict[coordId][tfID][speciesID]=[score]
		else:
			if(pval<cutOffVal):
				allTFDict[coordId][tfID][speciesID].append(score)
			#oldScore=allTFDict[coordId][tfID][speciesID][0]
			#if( oldScore<score):
				#allTFDict[coordId][tfID][speciesID]=[score, startPos, endPos, strand, pval, overlapDel]

#get the max difference
#its also possible to filter for a cutoff here, but we can do this downstream

outFile=open(outFileName, "w")
for coordKey in allTFDict:
    for tfKey in allTFDict[coordKey]:
        outArr=[]
        tempSpeciesDiffDict=dict()
        for speciesKey in allSpecies:
            #if the species is not in allSpecies, then liftover didn't map and we simply output NA
            if(speciesKey not in allTFDict[coordKey][tfKey]):
            #if(allTFDict[coordKey][tfKey][speciesKey]==[]):
                #outArr.extend(["NA"]*6)
                outArr.append("NA")
            else:
                #output avg score
                #if empty array, then output 0
                if(allTFDict[coordKey][tfKey][speciesKey]==[]):
                    avg_score=0
                else:
                    avg_score=np.mean(allTFDict[coordKey][tfKey][speciesKey])
                outArr.append(avg_score)
                tempSpeciesDiffDict[speciesKey]=avg_score
        
        #loop through all pairwise combinations, get difference
        for ind in range(len(pairwiseArr)):
            species1=pairwiseArr[ind][0]
            species2=pairwiseArr[ind][1]
            if((species1 not in tempSpeciesDiffDict) or (species2 not in tempSpeciesDiffDict)): 
                sumDiff="NA"
            else:
                sumDiff=tempSpeciesDiffDict[species2]-tempSpeciesDiffDict[species1]
            outArr.append(sumDiff)

        #write out the info for all species, one line per TF
        outFile.write(str(coordKey)+"\t"+str(tfKey)+"\t")
        for ind in range(len(outArr)):
            outFile.write(str(outArr[ind])+"\t")
        outFile.write("\n")


