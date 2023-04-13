#this script splits a combined fasta file for fimo for analyses

import sys

arg=sys.argv

fimoInputFileName=arg[1]
numCoordPerFile=int(arg[2])
outFileNameHeader=arg[3]

def writeFastaIterFile(fastaIterFileArr,outFileName):
    with open(outFileName, "w") as outFile:
        for line in fastaIterFileArr:
            fileID=line[0]
            sequence=line[1]
            outFile.write(fileID+"\n"+sequence+"\n")
        
#the fasta file must be sorted!

with open(fimoInputFileName, "r") as fimoInputFile:
    fileID=0
    fastaIterFileArr=[]
    prevCoordID="NA"
    numCoordSeen=-1

    for ind, l1 in enumerate(fimoInputFile):
        if(ind%2==0):     
            coordIdInfo=l1.split()[0]        
            coordIdInfoParsed=coordIdInfo.split("#")
            coordID=coordIdInfoParsed[0]
        else:

            if(coordID!=prevCoordID):
                numCoordSeen+=1
                prevCoordID=coordID

                #write out old array, reset if we have reached the number of unique coordinates per file         
                if(numCoordSeen==numCoordPerFile):
                    outFileName=outFileNameHeader+str(fileID)
                    writeFastaIterFile(fastaIterFileArr,outFileName)
                    #reset
                    numCoordSeen=0
                    fileID+=1
                    fastaIterFileArr=[]

            fastaSeq=l1.split()[0]
            fastaArrIter=[coordIdInfo, fastaSeq]
            fastaIterFileArr.append(fastaArrIter)

    #write out remaining lines at the end
    if(fastaIterFileArr!=[]):
        outFileName=outFileNameHeader+str(fileID)
        writeFastaIterFile(fastaIterFileArr,outFileName)

