#this script parses FIMO output and obtains the max difference between comparison files
#DON'T USE ${scriptdir}/fimo_get_max_diff_v2.py as there is a typo with calculating centerStartCoord and centerEndCoord, it didn't account for human deletion positions
#and the deletion overlap only had the situation of complete tf overlapping del
import sys

arg=sys.argv

fimoInputFileName=arg[1]
pairwiseCombFileName=arg[2]
outFileName=arg[3]
lookAtCenteredOnly=int(arg[4])

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
		
		
		#if( ((startPos<=centerStartCoord) and (centerEndCoord<=endPos)) or ((centerStartCoord<=startPos) and (endPos<=centerEndCoord)) or ((centerStartCoord<=startPos) and (centerEndCoord<=endPos)):
		if speciesID=="hg38":
			if( (endPos<=centerStartCoord) or (centerEndCoord<=startPos) ):
				overlapDel=0
			else:
				overlapDel=1
		else:
			if( (endPos<centerStartCoord) or (centerEndCoord<startPos) ):
				overlapDel=0
			else:
				overlapDel=1
		
		coordId=coordIdInfo[0]
		
		if(coordId not in allTFDict):
			allTFDict[coordId]=dict()
		if(tfID not in allTFDict[coordId]):
			allTFDict[coordId][tfID]=dict()
			allTFDict[coordId][tfID][speciesID]=[score, startPos, endPos, strand, pval, overlapDel]
		
		elif(speciesID not in allTFDict[coordId][tfID]):
			allTFDict[coordId][tfID][speciesID]=[score, startPos, endPos, strand, pval, overlapDel]
		else:
			oldScore=allTFDict[coordId][tfID][speciesID][0]
			if( oldScore<score):
				allTFDict[coordId][tfID][speciesID]=[score, startPos, endPos, strand, pval, overlapDel]

#get the max difference
#its also possible to filter for a cutoff here, but we can do this downstream
#first 2+6*3 columns are not pairwise info, next 3 are
outFile=open(outFileName, "w")
for coordKey in allTFDict:
    for tfKey in allTFDict[coordKey]:
        outArr=[]
        #sort to know output order
        for speciesKey in allSpecies:
            #if the species is not in allSpecies, then liftover didn't map and we simply output NA
            if(speciesKey not in allTFDict[coordKey][tfKey]):
                outArr.extend(["NA"]*6)
            else:
                outArr.extend(allTFDict[coordKey][tfKey][speciesKey])
 
        #loop through all pairwise combinations, get difference
        for ind in range(len(pairwiseArr)):
            species1=pairwiseArr[ind][0]
            species2=pairwiseArr[ind][1]
            if((species1 not in allTFDict[coordKey][tfKey]) or (species2 not in allTFDict[coordKey][tfKey])): 
                maxDiff="NA"
            else:
                maxDiff=allTFDict[coordKey][tfKey][species2][0]-allTFDict[coordKey][tfKey][species1][0]
            outArr.append(maxDiff)

        #write out the info for all species, one line per TF
        outFile.write(str(coordKey)+"\t"+str(tfKey)+"\t")
        for ind in range(len(outArr)):
            outFile.write(str(outArr[ind])+"\t")
        outFile.write("\n")


