#this script reads in the original coordinates and revises the coordinates to the chimp-human hybrid genome

import sys
arg = sys.argv


#new_del_arr[ind]=[human_gap_pos_chr, new_del_pos_start, new_del_pos_end, chimp_chr, chimp_chr_start_pos, chimp_chr_end_pos]

#new coordinates filename
chimpHybridCoordinatesFileName=arg[1]

#old coordinates
coordinatesToConvertFileName=arg[2]

#old coordinates with corresponding new coordinates
convertedCoordinatesOutFileName=arg[3]

human_gap_list=[]
coordIDConvertDict=dict()
with open(chimpHybridCoordinatesFileName, "r") as chimpHybridCoordinatesFile:
	for ind, l1 in enumerate(chimpHybridCoordinatesFile):
		#the last 3 coordinates will contain the original chimp coordinates
        #we use that as the id
		line=l1.split()		
		id=''.join(line[len(line)-3:len(line)])

		hybrid_pos_chr=line[0]
		hybrid_pos_start=int(line[1])
		hybrid_pos_end=int(line[2])
		
		coordIDConvertDict[id]=[hybrid_pos_chr, hybrid_pos_start, hybrid_pos_end]
		
#each line contains a the original chimp deletion/conserved coordinate, along with the conserved/deletion coordinates on the same line for conversion
newConvertedCoordArr=[]

with open(coordinatesToConvertFileName, "r") as coordinatesToConvertFile:
    for ind, l1 in enumerate(coordinatesToConvertFile):
        line=l1.split()
        print(ind)
        #the first 3 columns will contain the original chimp coordiantes
        orig_chimp_chr=line[0]
        orig_chimp_start=int(line[1])
        orig_chimp_end=int(line[2])
        id=''.join([orig_chimp_chr, str(orig_chimp_start), str(orig_chimp_end)])
        
        
        #obtain the chimp/human hybrid position for the overarching coordinate
        conversionCoord=coordIDConvertDict[id]
        conversionCoord_chr=conversionCoord[0]
        conversionCoord_start=conversionCoord[1]
        conversionCoord_end=conversionCoord[2]
        
        #indeces 6,7,8 contain the coordinates lying within the orig_chimp coordinates
        #coord_start is the original chimp coordinates
        coord_start=int(line[7])
        coord_end=int(line[8])
        coord_gap=coord_end-coord_start

        #get inversion information 
        invInfo=line[12]
        if(invInfo=="-"):
            coordOffset=orig_chimp_end-coord_end
            new_coord_start=coordOffset+conversionCoord_start
            new_coord_end=new_coord_start+coord_gap
        
        else:        
            coordOffset=coord_start-orig_chimp_start
        	
            #convert to chimp/human hybrid coordinates
            new_coord_start=coordOffset+conversionCoord_start
            new_coord_end=new_coord_start+coord_gap

        #store the coordinates
        newConvertedCoordArr.append(line)
        #insert the new converted coordinates
        #newConvertedCoordArr[len(newConvertedCoordArr)-1][9:9]=[orig_chimp_chr,new_coord_start, new_coord_end]
        #print(ind)
        #print([orig_chimp_chr,new_coord_start, new_coord_end])
        newConvertedCoordArr[len(newConvertedCoordArr)-1].extend([conversionCoord_chr, conversionCoord_start, conversionCoord_end, conversionCoord_chr,new_coord_start, new_coord_end])
        #print([orig_chimp_chr, conversionCoord_start, conversionCoord_end, orig_chimp_chr,new_coord_start, new_coord_end])
        #now the first four coordinates will be the "primary columns" - the columns that are the same between the deleted and conserved coordinates files

#write out the converted coordinates
#output will be the same as the original read file, except after the third coordinate (which contains the sequence overlapped)
#it will now contain the new chimp/human hybrid coordinate
with open(convertedCoordinatesOutFileName, "w") as convertedCoordinatesOutFile:
	for ind in range(len(newConvertedCoordArr)):
		for ind2 in range(len(newConvertedCoordArr[ind])):
			convertedCoordinatesOutFile.write(str(newConvertedCoordArr[ind][ind2])+"\t")
		convertedCoordinatesOutFile.write("\n")



