#this script parses FIMO output and obtains the max difference between comparison files
#v3 takes a file id encode TF map to match with FIMO
import sys

arg=sys.argv

fimoInputFileName=arg[1]
pValInd1=int(arg[2])
pValInd2=int(arg[3])
maxDiffInd=int(arg[4])
cutOffVal=float(arg[5])
useEncodeAnnotation=int(arg[6])
encodeFileName=arg[7]
encodeIDMapFileName=arg[8]
lookAtExpressedTFs=int(arg[9])
#this is a list of two files which contains the TF ids for species 1 and 2
expressedTFsFileName=arg[10]
FIMOTFIDMapFileName=arg[11]
outFileName=arg[12]

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
with open(fimoInputFileName, "r") as fimoInputFile:
	for ind, l1 in enumerate(fimoInputFile):
		l1=l1.split()
		coordId=l1[0]
		#convert TFId from FIMO into one thats comparable with encode
		TFIdInitial=l1[1]
		TFIdParsed=FIMOTFIDMap[TFIdInitial]
		
		maxDiff=float(l1[maxDiffInd])
		pVal1=float(l1[pValInd1])
		pVal2=float(l1[pValInd2])
		
		overlapDel1=int(l1[pValInd1+1])
		overlapDel2=int(l1[pValInd2+1])
		
		l1Info=l1[1:len(l1)]
		#the default is 0 for all columns
		if(coordId not in coordDict):
			coordDict[coordId]=[0]*len(l1Info)
		
		#we require at least one ref or alt to have a significant binding value AND that one or the other overlaps the actual deletion size 
		if( (pVal1<cutOffVal or pVal2<cutOffVal) and (overlapDel1==1 and overlapDel2==1) ):
			#if any element in the encode dict is not zero, then this means it overlaps with a TF
			#remember that a TF can have multiple files
			#look at maxDiffInd-1 due to looking at l1Info, which has the first coordinate removed
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
							if(abs(float(coordDict[coordId][maxDiffInd-1]))<abs(maxDiff)):
								coordDict[coordId]=l1Info
					else:
						if(abs(float(coordDict[coordId][maxDiffInd-1]))<abs(maxDiff)):
							coordDict[coordId]=l1Info
			
			#check if TF is expressed in both species
			elif(lookAtExpressedTFs==1):        
				if((TFIdInitial in expressedTFDictArr[0]) and (TFIdInitial in expressedTFDictArr[1])):
					if(abs(float(coordDict[coordId][maxDiffInd-1]))<abs(maxDiff)):
						coordDict[coordId]=l1Info        
			#if we don't want to use encode annotation, then simply keep the coordinate if its maxDiff is greater than the existing one  
			else:
				if(abs(float(coordDict[coordId][maxDiffInd-1]))<abs(maxDiff)):
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


