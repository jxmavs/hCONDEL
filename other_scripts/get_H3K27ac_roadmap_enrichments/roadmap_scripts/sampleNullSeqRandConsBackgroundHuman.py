#this script takes in a coordinate file and randomly samples deletion sites within the region
#unlike sampleNullSeqRandConsBackground.py, this version takes the corresponding human coordinates
#for this, we don't need to check for partial deletions, this is really the only modification

import random
import sys
import numpy as np

arg=sys.argv

coordinateFileName=arg[1]
maxTrys=int(arg[2])
consSizeDeviationPct=float(arg[3])
sampledOutCoordFileName=arg[4]

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
            sampledPosStart=random.randint(0, consLen-delLen)
            sampledPosEnd=sampledPosStart+delLen
            
            #if any of the sequence has been sampled already - i.e. annotated as a deletion, then can't sample the same region 
            if(coordInfoDict[consIDCoord][sampledPosStart:sampledPosEnd].any()):
                continue
            else:
                coordInfoDict[consIDCoord][sampledPosStart:sampledPosEnd]=True     
                sampledPosRelStart=sampledPosStart+consStartCoord
                sampledPosRelEnd=sampledPosEnd+consStartCoord
                allSampledCoord.append([consCoordChr, sampledPosRelStart, sampledPosRelEnd, idCoord]) 

                seqSampled=True 
                break  
            trys+=1 
        else:
            #this should rarely occur, but if it does, we have this, TODO: need to make a while loop here to keep sampling if say using consSizeDeviationPctIter fails also, keep extending until we get something
            print("maxTrys exceeded")    
 
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
