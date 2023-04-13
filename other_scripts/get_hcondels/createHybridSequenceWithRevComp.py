#this script creates the hybrid human/chimp sequence given the cleaned chimp/human hybrid coordinates
#output is fasta file containing hybrid sequence
#and a new metafile containing the coordinates on the hybrid sequence

import sys
from Bio.Seq import Seq
arg = sys.argv

#deletionCoordFile contains coordinates in format chimp pos, human pos
#and MUST be sorted by human coordinates, 

#deletionCoordFilePath contains the list of chimp coordinates to insert
deletionCoordFilePath=arg[1]
#this fasta file should match the coordinates in deletionCoordFilePath
deletionSequenceFilePath=arg[2]

#human sequence
humanFastaFilePath=arg[3]
humanChrLenFilePath=arg[4]

#out file names
hybrid_human_chimp_sequence_filename=arg[5]
hybrid_human_chimp_coordinates_metafilename=arg[6]

human_gap_list=[]


with open(deletionCoordFilePath, "r") as deletionCoordFile:
	for ind, l1 in enumerate(deletionCoordFile):
		line=l1.split()
		chimp_seq_chr=line[0]
		chimp_seq_chr_start=int(line[1])
		chimp_seq_chr_end=int(line[2])
		
		human_seq_chr=line[3]
		human_seq_chr_start=int(line[4])
		human_seq_chr_end=int(line[5])
		
		human_gap_list.append([chimp_seq_chr, chimp_seq_chr_start, chimp_seq_chr_end, human_seq_chr, human_seq_chr_start, human_seq_chr_end])

#read in chimp fasta sequences, this isn't the whole genome, but rather the fasta files from the deletion coordinates
chimp_fasta_seq=dict()

with open(deletionSequenceFilePath, "r") as deletionSequenceFile:
    for ind, l1 in enumerate(deletionSequenceFile):
        line=l1.split()
        lineEntry=line[0]
        #should be alternating between loading chr into dictionary and putting sequence into dictionary
        #each entry is a coordinate 
        if(">" == lineEntry[0]):
            chr_name=line[0][1:len(lineEntry)-2]
            chimp_fasta_seq[chr_name]=[]
            #get whehther or not sequence is inverted
            invInfo=line[0][len(lineEntry)-1]

        #if the sequence is inversion, then need to reverse complement the sequence
        else:
            if(invInfo=="-"):
                seq=Seq(lineEntry)
                revCompSeq=seq.reverse_complement()
                chimp_fasta_seq[chr_name]=list(revCompSeq)
            else:
                chimp_fasta_seq[chr_name]=list(lineEntry)
	
human_chr_len_dict=dict()
with open(humanChrLenFilePath, "r") as humanChrLenFile:
	for ind, l1 in enumerate(humanChrLenFile):
		line=l1.split()
		chr=line[0]
		chr_len=int(line[1])
		human_chr_len_dict[chr]=chr_len
		
human_chr_dict=dict()
with open(humanFastaFilePath, "r") as humanFastaFile:
	for ind, l1 in enumerate(humanFastaFile):
	
		
		line=l1.split()
		lineEntry=line[0]
		#this is the name for the chromosome
		if(">" == lineEntry[0]):
			#reset, get new chr
			chr_name=lineEntry[1:len(lineEntry)]			
			human_chr_dict[chr_name]=[0]*human_chr_len_dict[chr_name]
			currentPos=0
		else:
			#add on the sequence
			endPos=currentPos+len(lineEntry)
			human_chr_dict[chr_name][currentPos:endPos]=list(lineEntry)
			currentPos=endPos
        #should not happen
        if((""==lineEntry[0]) or (len(lineEntry)==0)):
            print(ind)

print("read hg38 file")			
		
		


#assemble the hybrid sequence, recalculate the base positions for the deletion coordinates
#each time you insert a sequence, the offset changes, offset is number you need to add/subtract to get correct position in the chimp/human hybrid genome
#IMPORTANT: need to make sure that human coordinates are sorted, otherwise nonsensical results will be outputted
#edge cases:  human deletion size extends into sequence skipped... should never occur?
offset=0
prev_pos_chr=human_gap_list[0][3]
new_chimp_arr=[0]*len(human_gap_list)

for ind in range(len(human_gap_list)):
	#print(ind)
    #initiation
	
	human_gap_pos_chr=human_gap_list[ind][3]
	

	#if we change chromosomes, then set offset to 0 again
	if(human_gap_pos_chr!=prev_pos_chr):
		offset=0
        prev_pos_chr=human_gap_pos_chr
		
	human_gap_pos_start=human_gap_list[ind][4]
	human_gap_pos_end=human_gap_list[ind][5]
	#length of sequence to skip in humans
	human_skipped_size=human_gap_pos_end-human_gap_pos_start
	
	#chimp position id
	
	chimp_chr=human_gap_list[ind][0]
	chimp_chr_start_pos=human_gap_list[ind][1]
	chimp_chr_end_pos=human_gap_list[ind][2]
	
	human_gap_size=chimp_chr_end_pos-chimp_chr_start_pos
	
	chimp_pos_id=chimp_chr+"|"+str(chimp_chr_start_pos)+"|"+str(chimp_chr_end_pos)
	
	human_insert_pos=human_gap_pos_start+offset
	
	#cut out/remove the skipped sequence
	if(human_skipped_size==0):
		#insert before specified position
		
		#get the chimp sequence
		chimp_seq=chimp_fasta_seq[chimp_pos_id]
	
		#insert the deleted chimp sequence
		human_chr_dict[human_gap_pos_chr][human_insert_pos:human_insert_pos]=chimp_seq
	
	else:
		#insert after specified position
		del human_chr_dict[human_gap_pos_chr][human_insert_pos:human_insert_pos+human_skipped_size]
		#get the chimp sequence
		chimp_seq=chimp_fasta_seq[chimp_pos_id]
	
		#insert the deleted chimp sequence
		human_chr_dict[human_gap_pos_chr][human_insert_pos:human_insert_pos]=chimp_seq
	
		
	#update new deletion position, keep track of meta information
	
	new_chimp_pos_start=human_insert_pos
	new_chimp_pos_end=human_insert_pos+len(chimp_seq)
	
	#store new chimp/human hybrid coordinates along with old chimp coordinates
	new_chimp_arr[ind]=[human_gap_pos_chr, new_chimp_pos_start, new_chimp_pos_end, chimp_chr, chimp_chr_start_pos, chimp_chr_end_pos]
	
	#update offset	
	offset=human_gap_size-human_skipped_size+offset
	

#human_chr_dict now contains the hybrid sequence

#write out the hybrid reference sequence
#follow same format as the other reference fasta sequences, 50 characters per line
with open(hybrid_human_chimp_sequence_filename, "w") as hybrid_chimp_sequence_file:
	
	chr_names_list=human_chr_dict.keys()
	for ind in range(len(chr_names_list)):
		chr=chr_names_list[ind]
		hybrid_chimp_sequence_file.write(">"+chr+"\n")
		
		#write out fastafile, each line contains 50 characters
		
		for ind2 in range(0, len(human_chr_dict[chr]),50):
			#join to combine all characters into single string
			partialseq=''.join(human_chr_dict[chr][ind2:ind2+50])
			hybrid_chimp_sequence_file.write(partialseq+"\n")

	
	
#write out the deletion coordinate positions/metadata
with open(hybrid_human_chimp_coordinates_metafilename, "w") as deletion_coordinates_metafile:
	for ind in range(len(new_chimp_arr)):
		for ind2 in range(len(new_chimp_arr[ind])):
			deletion_coordinates_metafile.write(str(new_chimp_arr[ind][ind2])+"\t")
		deletion_coordinates_metafile.write("\n")
	
	
			
			


