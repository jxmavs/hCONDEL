import sys

#this script takes a set of final filtered coordinates and creates the coordinates for the MPRA oligo design
#for partial deletions, I use the same deletion boundary as the deletions 

arg=sys.argv

coordinateFileName=arg[1]
chimpChrSizeFileName=arg[2]
humanChrSizeFileName=arg[3]
maxSeqLen=int(arg[4])
smallSeqOligoOutFileName=arg[5]
#large oligo file are the oligos which the deletion size is > maxSeqLen

smallSeqOligoOutFile=open(smallSeqOligoOutFileName,"w")

chimp_chr_sizes=dict()
with open(chimpChrSizeFileName, "r") as chimpChrSizeFile:
    for ind, l1 in enumerate(chimpChrSizeFile):
        l1=l1.split()
        chimp_chr_sizes[l1[0]]=int(l1[1])



human_chr_sizes=dict()
with open(humanChrSizeFileName, "r") as humanChrSizeFile:
    for ind, l1 in enumerate(humanChrSizeFile):
        l1=l1.split()
        human_chr_sizes[l1[0]]=int(l1[1])


#read in coordinates
with open(coordinateFileName, "r") as coordinateFile:
    for ind, l1 in enumerate(coordinateFile):
        l1=l1.split()
        chimpChr=l1[0]
        chimpSeqStartPos=int(l1[1])
        chimpSeqEndPos=int(l1[2])
        chimpSeqLen=chimpSeqEndPos-chimpSeqStartPos
        
        humanChr=l1[3]
        humanSeqStartPos=int(l1[4])
        humanSeqEndPos=int(l1[5])
        humanSeqLen=humanSeqEndPos-humanSeqStartPos
        #maxSeqLengthIter needs to account for human context, since there could be a gap on the human side as well
        minLengthIter=min(chimpSeqLen, humanSeqLen)

        #figure out how much to extend
        #if we have the deletion length greater than the max threshold
        #then simply center the oligo on the conserved position and output the 
        #oligo there, we don't output any human coordinates
        if(minLengthIter>maxSeqLen):
            
            
            newOligoCoordArr=[chimpChr, chimpSeqStartPos, chimpSeqEndPos, humanChr, humanSeqStartPos, humanSeqEndPos]
            
            for ind in range(len(newOligoCoordArr)):
                smallSeqOligoOutFile.write(str(newOligoCoordArr[ind])+"\t")

            #write out the rest of the input file
            for ind in range(len(l1)-1):
                smallSeqOligoOutFile.write(l1[ind]+"\t")
            smallSeqOligoOutFile.write(l1[len(l1)-1]+"\n")
 
            continue
    
        #otherwise, we extend the ends of the deletion for both human + chimp
        #we calculate rightExtendAmt and leftExtendAmt seperately since its possible that
        #we can't have exactly same bp of sequence on either end
        #but at least we can ensure that the amount of human/chimp bp or each respective end has the same bp
        rightExtendAmt=(maxSeqLen-minLengthIter)/2
        leftExtendAmt=(maxSeqLen-minLengthIter)-rightExtendAmt
   
        #set the beginning/end positions
        if(chimpSeqStartPos-leftExtendAmt<0):
            chimpOligoStartPos=0
        else:
            chimpOligoStartPos=chimpSeqStartPos-leftExtendAmt
        if(chimpSeqEndPos+rightExtendAmt>chimp_chr_sizes[chimpChr]):
            chimpOligoEndPos=chimp_chr_sizes[chimpChr]
        else:
            chimpOligoEndPos=chimpSeqEndPos+rightExtendAmt

        if(humanSeqStartPos-leftExtendAmt<0):
            humanOligoStartPos=0
        else:
            humanOligoStartPos=humanSeqStartPos-leftExtendAmt
        if(humanSeqEndPos+rightExtendAmt>human_chr_sizes[humanChr]):
            humanOligoEndPos=human_chr_sizes[humanChr]
        else:
            humanOligoEndPos=humanSeqEndPos+rightExtendAmt

       
        newOligoCoordArr=[chimpChr, chimpOligoStartPos, chimpOligoEndPos, humanChr, humanOligoStartPos, humanOligoEndPos]

        #write out the new oligo coordinates
        for ind in range(len(newOligoCoordArr)):
            smallSeqOligoOutFile.write(str(newOligoCoordArr[ind])+"\t")

        #write out the rest of the input file
        for ind in range(len(l1)-1):
            smallSeqOligoOutFile.write(l1[ind]+"\t")
        smallSeqOligoOutFile.write(l1[len(l1)-1]+"\n")

smallSeqOligoOutFile.close()
