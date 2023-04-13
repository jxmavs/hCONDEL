#this script parses FIMO output and obtains the max difference between comparison files
#v5 removes ENCODE and expressedTFs to speed up
import sys

arg=sys.argv

fimoInputFileName=arg[1]
pValInd1=int(arg[2])
pValInd2=int(arg[3])
maxDiffInd=int(arg[4])
cutOffVal=float(arg[5])
#this is a list of two files which contains the TF ids for species 1 and 2
FIMOTFIDMapFileName=arg[6]
outFileName=arg[7]


FIMOTFIDMap=dict()
with open(FIMOTFIDMapFileName, "r") as FIMOTFIDMapFile:
	for ind, l1 in enumerate(FIMOTFIDMapFile):
		l1=l1.split()
		TFIdInitial=l1[1]
		TFIdParsed=l1[3]
		FIMOTFIDMap[TFIdInitial]=TFIdParsed

#now it will be one TF entry per coordId
coordDict=dict()
with open(fimoInputFileName, "r") as fimoInputFile:
	for ind, l1 in enumerate(fimoInputFile):
		l1=l1.split()
		coordId=l1[0]
		#convert TFId from FIMO into one thats comparable with encode
		TFIdInitial=l1[1]
		TFIdParsed=FIMOTFIDMap[TFIdInitial]
		
		maxDiff=l1[maxDiffInd]
		pVal1=l1[pValInd1]
		pVal2=l1[pValInd2]
		
		overlapDel1=l1[pValInd1+1]
		overlapDel2=l1[pValInd2+1]
		
		l1Info=l1[1:len(l1)]
		#the default is 0 for all columns
		if(coordId not in coordDict):
			coordDict[coordId]=[0]*len(l1Info)
		
		#we require at least one ref or alt to have a significant binding value AND that one or the other overlaps the actual deletion size
		if(maxDiff!="NA"): 
			if( (float(pVal1)<cutOffVal or float(pVal2)<cutOffVal) and (int(overlapDel1)==1 and int(overlapDel2)==1) ):
				#remember that a TF can have multiple files
				#look at maxDiffInd-1 due to looking at l1Info, which has the first coordinate removed
				if(abs(float(coordDict[coordId][maxDiffInd-1]))<abs(float(maxDiff))):
					coordDict[coordId]=l1Info

#get the max difference
#its also possible to filter for a cutoff here, but we can do this downstream
#the ones without a maxdiff that passes the criteria will have a row of all zeros
with open(outFileName, "w") as outFile:
    for coordKey in coordDict:
        outFile.write(coordKey+"\t")
        for ind in range(len(coordDict[coordKey])):
            outFile.write(str(coordDict[coordKey][ind])+"\t")
        outFile.write("\n")


