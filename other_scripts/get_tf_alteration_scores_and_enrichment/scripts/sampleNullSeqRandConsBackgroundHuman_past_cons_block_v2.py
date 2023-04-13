#this script takes in a coordinate file and randomly samples deletion sites within the region
#unlike sampleNullSeqRandConsBackground.py, this version takes the corresponding human coordinates
#for this, we don't need to check for partial deletions, this is really the only modification

#made 10/25/21
#note that this doesn't accont for the situation where the deletion completely overlaps a conserved block - I don't have those

#10/28/21
#added keeping track of previously sampled bases in conserved sequences

import random
import sys
import numpy as np

arg=sys.argv

coordinateFileName=arg[1]
maxTrys=int(arg[2])
consSizeDeviationPct=float(arg[3])
sampledOutCoordFileName=arg[4]
chrSizeFileName=arg[5]
prevDrawnInfoFileName=arg[6]

chr_sizes=dict()
with open(chrSizeFileName, "r") as chrSizeFile:
	for ind, l1 in enumerate(chrSizeFile):
		#print(l1[0])
		#print(l1[1])
		l1=l1.split()
		chr_sizes[l1[0]]=int(l1[1])


coordInfoDict=dict()
allIDCoords=[]
with open(coordinateFileName, "r") as coordinateFile:
    for ind, l1 in enumerate(coordinateFile):
        l1=l1.split()

        consCoordChr=l1[0]
        consStartCoord=int(l1[1])
        consEndCoord=int(l1[2])

        origConsStartCoord=int(l1[4])
        origConsEndCoord=int(l1[5])
        #randomly sample 
        delStartCoord=int(l1[7])
        delEndCoord=int(l1[8])
   
        consLen=consEndCoord-consStartCoord
        
        #for partial deletions, we only consider the deletion length which overlies the entire conserved region
        #if((origConsStartCoord>=delStartCoord) and (origConsEndCoord>=delEndCoord)): 
        #    delLen=delEndCoord-origConsStartCoord
        #elif((origConsStartCoord<=delStartCoord) and (origConsEndCoord<=delEndCoord)):
        #    delLen=origConsEndCoord-delStartCoord
        #else:
        #    delLen=delEndCoord-delStartCoord
        delLen=delEndCoord-delStartCoord
        #if(delLen>consLen):
        #    delLen=consLen

        idCoord=l1[9]
      
        consIDCoord=''.join([consCoordChr, str(consStartCoord), str(consEndCoord)]) 
        if(consIDCoord not in coordInfoDict):
            coordInfoDict[consIDCoord]=np.array([False]*consLen)
        
        allIDCoords.append((consCoordChr, consStartCoord, consEndCoord, consLen, delLen, idCoord))

allCoordData=np.array(allIDCoords, dtype=[('consCoordChr', 'a100'), ('consStartCoord','i4'), ('consEndCoord', 'i4'), ('consLen', 'i4'), ('delLen', 'i4'), ('idCoord', 'a1000')  ])

#if sequence is previously drawn, then prevent from sampling again 
if prevDrawnInfoFileName!="NA":
	with open(prevDrawnInfoFileName, "r") as prevDrawnInfoFile:
		for ind, l1 in enumerate(prevDrawnInfoFile):	
						
			l1=l1.split()
			consCoordChr=l1[0]
			consStartCoord=int(l1[1])
			consEndCoord=int(l1[2])
			
			origConsStartCoord=int(l1[4])
			origConsEndCoord=int(l1[5])
			
			sampledDelStartCoord=int(l1[7])
			sampledDelEndCoord=int(l1[8])
			
			consLen=consEndCoord-consStartCoord
			
			idCoord=l1[9]
			
			initialSampledPosStart=sampledDelStartCoord-consStartCoord
			initialSampledPosEnd=sampledDelEndCoord-consStartCoord
			
			#check if you go over conserved block area
			if initialSampledPosEnd>consLen:
				sampledPosEnd=consLen
				sampledPosStart=initialSampledPosStart
			elif initialSampledPosStart<0:
				sampledPosEnd=initialSampledPosEnd
				sampledPosStart=0
			else:
				sampledPosEnd=initialSampledPosEnd
				sampledPosStart=initialSampledPosStart
						
			
			consIDCoord=''.join([consCoordChr, str(consStartCoord), str(consEndCoord)])
			if(consIDCoord in coordInfoDict):
				coordInfoDict[consIDCoord][sampledPosStart:sampledPosEnd]=True



allSampledCoord=[]

for ind in range(len(allCoordData)):
	#print(ind)
	consCoordChr=allCoordData[ind]['consCoordChr']
	consStartCoord=allCoordData[ind]['consStartCoord']
	consEndCoord=allCoordData[ind]['consEndCoord']
	consIDCoord=''.join([consCoordChr, str(consStartCoord), str(consEndCoord)])
	
	idCoord=allCoordData[ind]['idCoord']
	delLen=allCoordData[ind]['delLen']
	consLen=allCoordData[ind]['consLen']
	
	#origConsLen stores the consLen in case we need to draw another conserved sequence to sample the deletion site
	origConsLen=consLen
	consSizeDeviationPctIter=consSizeDeviationPct
	#keep sampling until you can get a position, 
	#it could be impossible to get a position, if we exceed maxTrys, sample the deletion in another place
	#that preferentially matches the size/attributes of the original position
	seqSampled=False
	while not(seqSampled):
		    
		trys=1
		while trys<=maxTrys:    
			
			#remember that sampled sequence is bed, assume that delLen is always greater than 1
			if np.random.randint(0, 2)==0:
				initialSampledPosStart=random.randint(0-(delLen-1), consLen-1)
				initialSampledPosEnd=initialSampledPosStart+delLen
			else:
				initialSampledPosEnd=random.randint(1, consLen+(delLen-1))
				initialSampledPosStart=initialSampledPosEnd-delLen
			
			#check if it goes past chromosome, resample
			if (initialSampledPosEnd+consStartCoord>chr_sizes[consCoordChr]) or (initialSampledPosStart+consStartCoord<0):
				trys+=1
				continue
			#check if you go over conserved block area
			if initialSampledPosEnd>consLen:
				sampledPosEnd=consLen
				sampledPosStart=initialSampledPosStart		
			elif initialSampledPosStart<0:
				sampledPosEnd=initialSampledPosEnd
				sampledPosStart=0	
			else:
				sampledPosEnd=initialSampledPosEnd
				sampledPosStart=initialSampledPosStart
			
			#if any of the sequence has been sampled already - i.e. annotated as a deletion, then can't sample the same region 
			if(coordInfoDict[consIDCoord][sampledPosStart:sampledPosEnd].any()):
				trys+=1
				continue
			else:
				coordInfoDict[consIDCoord][sampledPosStart:sampledPosEnd]=True     
				sampledPosRelStart=initialSampledPosStart+consStartCoord
				sampledPosRelEnd=initialSampledPosEnd+consStartCoord
				allSampledCoord.append([consCoordChr, sampledPosRelStart, sampledPosRelEnd, idCoord]) 
				
				seqSampled=True 
				break  
		else:
			#this should rarely occur, but if it does, we have this, TODO: need to make a while loop here to keep sampling if say using consSizeDeviationPctIter fails also, keep extending until we get something, changed 10/29/21 to just exit
			sys.exit("Max Trys Exceeded")    
			
			#we will draw another seqeunce with similar sequence properties as the sequence attempted to be sampled
			seqLowBound=origConsLen-consSizeDeviationPctIter*(origConsLen)
			seqUpBound=origConsLen+consSizeDeviationPctIter*(origConsLen)
			
			subSetInd=( (allCoordData['consLen']>=seqLowBound) & (allCoordData['consLen']<=seqUpBound) ).nonzero()[0]
			randInd=np.random.randint(0, len(subsetInd)-1)
			sampleInd=subsetInd[randInd]
			
			#note how we change every parameter EXCEPT delLen 
			consCoordChr=allCoordData[sampleInd]['consCoordChr']
			consStartCoord=allCoordData[sampleInd]['consStartCoord']
			consEndCoord=allCoordData[sampleInd]['consEndCoord']
			consIDCoord=''.join([consCoordChr, str(consStartCoord), str(consEndCoord)])
			
			idCoord=allCoordData[sampleInd]['idCoord']
			consLen=allCoordData[sampleInd]['consLen']
			
			
			#if we fail to sample again, then next time we'll have a larger consSizeDeviation to work with               
			consSizeDeviationPctIter=consSizeDeviationPctIter+0.01
			#generate another idCoord to sample
			#randIDCoordInd=random.randomint()      


with open(sampledOutCoordFileName, "w") as sampledOutCoordFile: 
    for ind in range(len(allSampledCoord)):
        for ind2 in range(len(allSampledCoord[ind])):
            sampledOutCoordFile.write(str(allSampledCoord[ind][ind2])+"\t")
        sampledOutCoordFile.write("\n")
