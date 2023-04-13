#this script takes an output file from bedtools closest and formats it to see if the distances are less than certain thresholds

import sys

arg = sys.argv

bedtoolsClosestDistOutputFileName=arg[1]
distInterestFilterFileName=arg[2]
outFileName=arg[3]

outFile=open(outFileName, "w")

distInterestFilterArr=[]
#read in distFilters
with open(distInterestFilterFileName, "r") as distInterestFilterFile:
	for ind, l1 in enumerate(distInterestFilterFile):
		l1=l1.split()
		distFilter=float(l1[0])
		distInterestFilterArr.append(distFilter)

#read in coordinates
with open(bedtoolsClosestDistOutputFileName, "r") as bedtoolsClosestDistOutputFile:
	for ind, l1 in enumerate(bedtoolsClosestDistOutputFile):					 
		l1=l1.split()
		seqId=l1[0]
		closestDist=int(l1[1])
		
		outFile.write(seqId+"\t")
		for distFilterInd in range(len(distInterestFilterArr)):
			distFilter=distInterestFilterArr[distFilterInd]
			if closestDist!=-1 and closestDist<distFilter:
				outFile.write("1\t")
			else:
				outFile.write("0\t")
		
		outFile.write(str(closestDist)+"\n")	

outFile.close()

 
