import sys

#this script takes a set of coordinates and extends them by ${extendBP}, if it encounters a 
#less than start site, or goes over the actual chromosome size, it will readjust
arg = sys.argv

coordinateFileName = arg[1]
chrSizeFileName=arg[2]
extendBP_up=int(arg[3])
extendBP_down=int(arg[4])
outFileName=arg[5]

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
		chr=l1[0]
		startPos=int(l1[1])
		endPos=int(l1[2])
		strand=l1[3]
		
		if((strand=="1") or (strand=="+")):
				
			if(startPos-extendBP_up<0):
				newStartPos=0
			else:
				newStartPos=startPos-extendBP_up
							
			if(endPos+extendBP_down>chr_sizes[chr]):
				newEndPos=chr_sizes[chr]
			else:
				newEndPos=endPos+extendBP_down
		else:
			
			if(startPos-extendBP_down<0):
				newStartPos=0
			else:
				newStartPos=startPos-extendBP_down
			
			if(endPos+extendBP_up>chr_sizes[chr]):
				newEndPos=chr_sizes[chr]
			else:
				newEndPos=endPos+extendBP_up	
	
		outFile.write(chr+"\t"+str(newStartPos)+"\t"+str(newEndPos)+"\t")
		for ind in range(len(l1)-1):
			outFile.write(l1[ind]+"\t")
		outFile.write(l1[len(l1)-1]+"\n")

outFile.close()
