import sys

#this is a general script that extends coordinates centered on the site of interest
#and outputs extended coordinates as well as the original deletion coordinates
#${scriptdir}/extendCoordinatesCHIPV3.py

arg = sys.argv

coordinateFileName=arg[1]
chrSizeFileName=arg[2]
maxBPWindow=int(arg[3])
outFileName=arg[4]

outFile=open(outFileName,"w")

#read in chromosome sizes

chr_sizes=dict()
with open(chrSizeFileName, "r") as chrSizeFile:
	for ind, l1 in enumerate(chrSizeFile):
		#print(l1[0])
		#print(l1[1])
		l1=l1.split()
		chr_sizes[l1[0]]=int(l1[1])


#read in coordinates
with open(coordinateFileName, "r") as coordinateFile:
    for ind, l1 in enumerate(coordinateFile):
        l1=l1.split()

        chrName=l1[0]
        chrStartPos=int(l1[1])
        chrEndPos=int(l1[2])
        invInfo=l1[3]
        idCoord=l1[4]
        
        leftExtendAmt=maxBPWindow/2
        rightExtendAmt=maxBPWindow-leftExtendAmt
        
        #chrMiddleStartPos=(chrEndPos-chrStartPos)/2+chrStartPos
        #take into account inversions in extending
        if(invInfo=="+"):
            #set the beginning/end positions
            if(chrStartPos-leftExtendAmt<0):
                chrExtendStartPos=0
            else:
                chrExtendStartPos=chrStartPos-leftExtendAmt
            
            if(chrEndPos+rightExtendAmt>chr_sizes[chrName]):
                chrExtendEndPos=chr_sizes[chrName]
            else:
                chrExtendEndPos=chrEndPos+rightExtendAmt
            
        else:
            #set the beginning/end positions
            if(chrStartPos-rightExtendAmt<0):
                chrExtendStartPos=0
            else:
                chrExtendStartPos=chrStartPos-rightExtendAmt
            
            if(chrEndPos+leftExtendAmt>chr_sizes[chrName]):
                chrExtendEndPos=chr_sizes[chrName]
            else:
                chrExtendEndPos=chrEndPos+leftExtendAmt

         
        '''       
        
        if(chrEndPos-chrStartPos>maxBPWindow):
            chrExtendStartPos=chrStartPos
            chrExtendEndPos=chrEndPos
        #elif(maxBPWindow!=0):      
        else:
            leftExtendAmt=maxBPWindow/2
            rightExtendAmt=maxBPWindow-leftExtendAmt

            chrMiddleStartPos=(chrEndPos-chrStartPos)/2+chrStartPos
            #set the beginning/end positions
            if(chrMiddleStartPos-leftExtendAmt<0):
                chrExtendStartPos=0
            else:
                chrExtendStartPos=chrMiddleStartPos-leftExtendAmt

            if(chrMiddleStartPos+rightExtendAmt>chr_sizes[chrName]):
                chrExtendEndPos=chr_sizes[chrName]
            else:
                chrExtendEndPos=chrMiddleStartPos+rightExtendAmt
        '''
        '''
        leftExtendAmt=maxBPWindow/2
        rightExtendAmt=maxBPWindow-leftExtendAmt

        chrMiddleStartPos=(chrEndPos-chrStartPos)/2+chrStartPos


        #set the beginning/end positions
        if(chrMiddleStartPos-leftExtendAmt<0):
            chrExtendStartPos=0
        else:
            chrExtendStartPos=chrMiddleStartPos-leftExtendAmt

        if(chrMiddleStartPos+rightExtendAmt>chr_sizes[chrName]):
            chrExtendEndPos=chr_sizes[chrName]
        else:
            chrExtendEndPos=chrMiddleStartPos+rightExtendAmt
        '''

        outFile.write(chrName+"\t"+str(chrExtendStartPos)+"\t"+str(chrExtendEndPos)+"\t"+chrName+"\t"+str(chrStartPos)+"\t"+str(chrEndPos)+"\t"+idCoord+"\n")


outFile.close()
