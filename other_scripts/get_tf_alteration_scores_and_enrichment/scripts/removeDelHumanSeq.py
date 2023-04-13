#remove human annotated deleted sequence from fasta

import sys

arg=sys.argv

fastaFileName=arg[1]
outFileName=arg[2]

outFile=open(outFileName, "w")

with open(fastaFileName, "r") as fastaFile:
	for ind, l1 in enumerate(fastaFile):
		l1=l1.split()[0]
		
		if ind%2==0:
			coordIdInfo=l1
			
			coordIdInfo=coordIdInfo.split("#")
			speciesInfo=coordIdInfo[1]
			speciesInfo=speciesInfo.split("|")
			
			#first 3 coordinates are the overarching sequence, next 3 are the central coordinates, remember they are bed coordinates, so centerEndPos does not include a base
			centerStartCoord=int(speciesInfo[5])-int(speciesInfo[2])
			centerEndCoord=int(speciesInfo[6])-int(speciesInfo[2])
			
			#revise coordinate id
			#deletion end position is now the same as start, for FIMO
			delLen=int(speciesInfo[6])-int(speciesInfo[5])
			
			speciesInfo[5]=speciesInfo[6]
			speciesInfo[2]=str(int(speciesInfo[2])+delLen)
			
			speciesInfoNew="|".join(speciesInfo)
			coordIdInfo[1]=speciesInfoNew
			coordIdInfoNew="#".join(coordIdInfo)
			
		else:
			sequence=l1	
			#remove deleted sequence
			revised_sequence=sequence[0:centerStartCoord]+sequence[centerEndCoord:len(sequence)]
			
			outFile.write(coordIdInfoNew+"\n"+revised_sequence+"\n")

outFile.close()

