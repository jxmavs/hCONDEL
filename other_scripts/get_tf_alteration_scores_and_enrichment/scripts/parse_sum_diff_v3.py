#this script parses FIMO output and obtains the sum difference between comparison files
#v2 takes encodeIDMapFileName and maps fileID to TFName
import sys

arg=sys.argv

fimoInputFileName=arg[1]
sumDiffInd=int(arg[2])
useEncodeAnnotation=int(arg[3])
encodeFileName=arg[4]
encodeIDMapFileName=arg[5]
lookAtExpressedTFs=int(arg[6])
#this is a list of two files which contains the TF ids for species 1 and 2
expressedTFsFileName=arg[7]
FIMOTFIDMapFileName=arg[8]
outFileName=arg[9]

#file id, gene name
fileIDTFDict=dict()
with open(encodeIDMapFileName, "r") as encodeIDMapFile:
    for ind, l1 in enumerate(encodeIDMapFile):
        l1=l1.split()
        fileID=l1[0]
        TFName=l1[1]
        fileIDTFDict[fileID]=TFName

FIMOTFIDMap=dict()
with open(FIMOTFIDMapFileName, "r") as FIMOTFIDMapFile:
    for ind, l1 in enumerate(FIMOTFIDMapFile):
        l1=l1.split()
        TFIdInitial=l1[1]
        TFIdParsed=l1[3]
        FIMOTFIDMap[TFIdInitial]=TFIdParsed

#encodeArr contains the files which lists whether or not the coordinate lies in a peak region
#encodeFile is a matrix with a value greater than 1 if there is an intersection
encodeDict=dict()
with open(encodeFileName, "r") as encodeFile:
    
    TFHeaderLine=encodeFile.readline()
    TFInfo=TFHeaderLine.split() 
    allTFNames=[]

    #remove the fileID from the TF names
    for ind in range(len(TFInfo)):
        fileID=TFInfo[ind].split("_")[1]
        TFName=fileIDTFDict[fileID]
        allTFNames.append(TFName)

    #now, look at the rest of the lines    
    for ind, l1 in enumerate(encodeFile):
        l1=l1.split()
        coordId=l1[0]
        encodeDict[coordId]=dict()

        TFMarkInfo=l1[1:len(l1)]
        for ind2 in range(len(allTFNames)):
            TFId=allTFNames[ind2]

            if(TFId not in encodeDict):
                encodeDict[coordId][TFId]=[]
            #append on info since we can have multiple experiments with the same TF, we look at all of them
            encodeDict[coordId][TFId].append(int(TFMarkInfo[ind2]))

#this is just a list of TFs expressed in the cell type
expressedTFDictArr=[dict(), dict(), dict()]
#expressedTFsFileName contains a a list of files, which contain the TFs expressed for that particular species - one line for species 1, another for species 2

with open(expressedTFsFileName, "r") as expressedTFsFile:
    ind=0
    for line in expressedTFsFile:
        speciesTFFileName=line.split()[0]

        with open(speciesTFFileName, "r") as speciesTFFile:

            for line2 in speciesTFFile:
                
                TFId=line2.split()[0]
                expressedTFDictArr[ind][TFId]=1
        ind+=1
#now it will be one TF entry per coordId
coordDict=dict()
#mark if the coordinate has a TF annotation
coordTFAnnoteDict=dict()
with open(fimoInputFileName, "r") as fimoInputFile:
	for ind, l1 in enumerate(fimoInputFile):
		l1=l1.split()
		coordId=l1[0]
		#print(l1) 
		TFIdInitial=l1[1]
		TFIdParsed=FIMOTFIDMap[TFIdInitial]
		
		if coordId not in coordTFAnnoteDict:
			coordTFAnnoteDict[coordId]=0
		
		sumDiff=l1[sumDiffInd]
		#set to 0, so it won't be set to any value in the if statements below, but this should never be NA unless I'm specifically focusing on macaque/human macaque/chimp differences..
		if(sumDiff=="NA"):
			sumDiff=0
		else:
			sumDiff=float(sumDiff)
			coordTFAnnoteDict[coordId]=1
		
		l1Info=l1[1:len(l1)]
		#the default is 0 for all columns
		if(coordId not in coordDict):
			coordDict[coordId]=[0]*len(l1Info)
		
		#look at sumDiffInd-1 due to us looking at l1Info, which does not contain the first element 
		if(useEncodeAnnotation==1):
			#if complex, require both of the TFs to have peaks
			if("::" in TFIdParsed):
				TFIDParsedArr=TFIdParsed.split("::")
				for ind in range(len(TFIDParsedArr)):
					TFIdParsed=TFIDParsedArr[ind]
					if(TFIdParsed in encodeDict[coordId]):
						encodePeakFound=not(all([ iterItem == 0 for iterItem in encodeDict[coordId][TFIdParsed] ]) )
						if not(encodePeakFound):
							break
					else:
						encodePeakFound=False
						break
			else:
				if(TFIdParsed in encodeDict[coordId]):
					encodePeakFound=not(all([ iterItem == 0 for iterItem in encodeDict[coordId][TFIdParsed] ]) )
				else:
					encodePeakFound=False
						
			if encodePeakFound:
				if(lookAtExpressedTFs==1):
					if((TFIdInitial in expressedTFDictArr[0]) and (TFIdInitial in expressedTFDictArr[1])):
						if(abs(float(coordDict[coordId][sumDiffInd-1]))<abs(sumDiff)):
							coordDict[coordId]=l1Info
				else:
					if(abs(float(coordDict[coordId][sumDiffInd-1]))<abs(sumDiff)):
						coordDict[coordId]=l1Info
		
		#check if TF is expressed in both species
		elif(lookAtExpressedTFs==1):        
			if((TFIdInitial in expressedTFDictArr[0]) and (TFIdInitial in expressedTFDictArr[1])):
				if(abs(float(coordDict[coordId][sumDiffInd-1]))<abs(sumDiff)):
					coordDict[coordId]=l1Info        
		#if we don't want to use encode annotation, then simply keep the coordinate if its sumDiff is greater than the existing one  
		else:
			if(abs(float(coordDict[coordId][sumDiffInd-1]))<abs(sumDiff)):
				coordDict[coordId]=l1Info

#get the sum difference
#its also possible to filter for a cutoff here, but we can do this downstream
#the ones without a sumdiff that passes the criteria will have a row of all zeros
with open(outFileName, "w") as outFile:
    for coordKey in coordDict:
        outFile.write(coordKey+"\t")
        for ind in range(len(coordDict[coordKey])):
            outFile.write(str(coordDict[coordKey][ind])+"\t")
        outFile.write(str(coordTFAnnoteDict[coordId]))
        outFile.write("\n")


