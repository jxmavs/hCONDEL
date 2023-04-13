#this script parses FIMO output and obtains the max difference between comparison files

import sys

arg=sys.argv

memeFileName=arg[1]
outFileHeader=arg[2]

header_lines=[]

headerFound=False
outFileIter="not opened"
with open(memeFileName, "r") as memeFile:
	for ind, l1 in enumerate(memeFile):
		#keep the first few lines before 
		if (not headerFound):
			if "MOTIF" in l1:
				headerFound=True
			else:
				header_lines.append(l1)
				#ignore lines below
				continue
		
		if "MOTIF" in l1:
			if outFileIter!="not opened":
				outFileIter.close()
			
			l1_split=l1.split()
			TF_name=l1_split[2]
			outFileIter=open(outFileHeader+"."+TF_name+".txt", "w")
			
			#write out header lines first
			for ind in range(len(header_lines)):
				outFileIter.write(header_lines[ind])
			#write out motif header
			outFileIter.write(l1)	
		else:
			#write out rest of file
			outFileIter.write(l1)

#close last meme file
outFileIter.close()

