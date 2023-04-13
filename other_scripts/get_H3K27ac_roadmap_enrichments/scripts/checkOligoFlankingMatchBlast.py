import sys

#this script takes a raw sam aligned file and checks to see the number of mismatches/gaps for readName==refName

arg = sys.argv

inputFileName=arg[1]
outFileName=arg[2]

totalNumMismatch=0
totalNumGapOpen=0
totalNumGaps=0

with open(inputFileName, "r") as inputFile:
    for ind, l1 in enumerate(inputFile):
        l1=l1.split()
        readName=l1[0]
        qLen=int(l1[3])
        startPos=int(l1[5])
        endPos=int(l1[6])

        numMismatch=int(l1[9])
        numGapOpen=int(l1[10])
        numGaps=int(l1[11])
               
        if ind==0:
            qAlignArr=[0]*qLen
            
        #if an alignment is split into multiple alignments, then only keep the alignment that already filled the overlapping positions, as the blast file is sorted from best to worse alignment
        if 1 in qAlignArr[(startPos-1):endPos]:
            continue   
        else:
  
            totalNumMismatch+=numMismatch
            totalNumGapOpen+=numGapOpen
            totalNumGaps+=numGaps
         
            qAlignArr[(startPos-1):endPos]=[1]*(endPos-startPos+1)
            
numUnalignedBp=0
for ind in range(len(qAlignArr)):
    if qAlignArr[ind]==0:
        numUnalignedBp+=1

with open(outFileName, "w") as outFile:
    outFile.write(readName+"\t"+str(totalNumMismatch)+"\t"+str(totalNumGapOpen)+"\t"+str(totalNumGaps)+"\t"+str(numUnalignedBp)+"\t"+str(numUnalignedBp+totalNumMismatch+totalNumGaps)+"\n")
