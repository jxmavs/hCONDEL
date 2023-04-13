import sys
import numpy as np
import random

random.seed(10)

arg=sys.argv

closestDistFileName=arg[1]
outFileName=arg[2]

closestDistOutputDict=dict()
with open(closestDistFileName, "r") as closestDistFile:
	for ind, l1 in enumerate(closestDistFile):
		l1=l1.split()
		seqName=l1[4]
		tpm=float(l1[9])
		closestDistOutputDict[seqName]=l1
		if seqName in closestDistOutputDict:
			tpm_prev=float(closestDistOutputDict[seqName][9])
			if tpm>tpm_prev:
				closestDistOutputDict[seqName]=l1
			elif tpm==tpm_prev:
				if random.randint(0, 1)==0:
					closestDistOutputDict[seqName]=l1
			else:
				a=1


outFile=open(outFileName, "w")
for seqName in closestDistOutputDict:
	l1=closestDistOutputDict[seqName]
	for ind in range(len(l1)-1):
		outFile.write(l1[ind]+"\t")
	outFile.write(l1[len(l1)-1]+"\n")

outFile.close()
