#trim sequence to optmize for FIMO run

import sys

arg=sys.argv

fastaFileName=arg[1]
motifLength=int(arg[2])
outFileName=arg[3]

#if want buffer at the end of sequence, not required
bufferLength=0

outFile=open(outFileName, "w")

with open(fastaFileName, "r") as fastaFile:
	for ind, l1 in enumerate(fastaFile):
		l1=l1.split()[0]
		
		if ind%2==0:		
			coordIdInfo=l1
			coordIdInfo=coordIdInfo.split("#")
			speciesInfo=coordIdInfo[1]
			speciesInfo=speciesInfo.split("|")
			
			seqBoundaryStartPos=int(speciesInfo[2])
			
			#first 3 coordinates are the overarching sequence, next 3 are the central coordinates, remember they are bed coordinates, so centerEndPos does not include a base
			centerStartCoord=int(speciesInfo[5])-seqBoundaryStartPos
			centerEndCoord=int(speciesInfo[6])-seqBoundaryStartPos
			
			if(centerStartCoord-(motifLength-1)-bufferLength<0):
				trimStartSeq=0
			else:
				trimStartSeq=centerStartCoord-(motifLength-1)-bufferLength
			
			if(centerEndCoord+(motifLength-1)+bufferLength>(int(speciesInfo[3])-seqBoundaryStartPos)):
				trimEndSeq=int(speciesInfo[3])-seqBoundaryStartPos
			else:
				trimEndSeq=centerEndCoord+(motifLength-1)+bufferLength
			
			#revise coordinate id
			#change boundaries of overlapping sequence
			speciesInfo[2]=str(trimStartSeq+seqBoundaryStartPos)
			speciesInfo[3]=str(trimEndSeq+seqBoundaryStartPos)
			
			speciesInfoNew="|".join(speciesInfo)
			coordIdInfo[1]=speciesInfoNew
			coordIdInfoNew="#".join(coordIdInfo)
			
		else:
			sequence=l1
			#trim sequence
			revised_sequence=sequence[trimStartSeq:trimEndSeq]
			
			outFile.write(coordIdInfoNew+"\n"+revised_sequence+"\n")

outFile.close()

