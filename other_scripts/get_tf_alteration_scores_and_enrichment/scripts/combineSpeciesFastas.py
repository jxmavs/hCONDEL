#this script aggregates files and combines it, writes it out together

import sys

arg=sys.argv

inputFastaListFileName=arg[1]
outFileName=arg[2]

outFile=open(outFileName, "w")

aggregated_files=[]
with open(inputFastaListFileName, "r") as inputFastaListFile:
	for ind, l1 in enumerate(inputFastaListFile):
		inputFastaListFileIterName=l1.split()[0]
		aggregated_files.append([])
		
		with open(inputFastaListFileIterName, "r") as inputFastaListFileIter:		
			for ind2, l2 in enumerate(inputFastaListFileIter):
				file_iter_line=l2.split()[0]
				aggregated_files[ind].append(file_iter_line)

#print(aggregated_files)
#combine fastas, write out files
for ind2 in range(0, len(aggregated_files[0]), 2):
	for ind in range(len(aggregated_files)-1):
		outFile.write(aggregated_files[ind][ind2]+"\n")
		outFile.write(aggregated_files[ind][ind2+1]+"\n")
	outFile.write(aggregated_files[len(aggregated_files)-1][ind2]+"\n")
	outFile.write(aggregated_files[len(aggregated_files)-1][ind2+1]+"\n")
outFile.close()

