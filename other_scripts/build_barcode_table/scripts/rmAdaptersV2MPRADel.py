## Dustin Griesemer
##
## 02/15/15
##
## rmAdapters.py
##
## usage: python rmAdapters.py fastq_file_path -o output_file_path -d discard_file_path [-Lbarcode Lbarcode] [-Rbarcode Rbarcode] [-Ladapter Ladapter] [-seqxN] [-barxN]
##
## example: python rmAdapters.py /cluster_path/sequencing/AALE9_2/fastq/AALE9.AACTTGAC.extendedFrags.fastq -o trimmed.expos.fastq -d discarded.expos.fastq 
##
## Removes adapters from an MPRA fastq file, cataloging barcodes in the process. Outputs an updated fastq file, with the adapters trimmed and the barcode in the sequence identifier. If the -seqxN option is used, sequences with N's are removed. If the -barxN option is used, sequences with N's in the barcode position are removed.

#adapted from Dustin's rmAdapters code
import sys
import os
import argparse
import fuzzysearch as fs
from Bio import SeqIO
from Bio.Seq import Seq
import Levenshtein as lv
import gzip

#if we don't want to use fuzzysearch... function still needs to be fixed
def Radapter_trim(seq,Lbar,Rbar,Nbar=20,Lbar_dist=3,Rbar_dist=3,Bar_clamp=2):
#Takes in sequence in string form, outputs index at which trimming should occur and barcode (or error) sequence
#Avoids using fuzzysearch for efficiency
#Searches for an exact match to 10bp on either the left or right side of the barcode.
#Finds Lbar & Rbar based on position and performs Levenshtein distance check.
#Nbar: expected length of barcode
#Finds Barcode, requiring an exact match of Bar_clamp bases on either side.
	if Lbar[-10:] in seq:
		Lbar_pos=seq.find(Lbar[-10:])-len(Lbar)+10
		Bar_pos=Lbar_pos+len(Lbar)
		Rbar_pos=Bar_pos+Nbar
	elif Rbar[0:10] in seq:
		Rbar_pos=seq.find(Rbar[0:10])
		Bar_pos=Rbar_pos-Nbar
		Lbar_pos=Bar_pos-len(Lbar)
	else:
		return 'Lbar/Rbar_not_found',0 #Couldn't find an exact match on either left or right of barcode
	Lbar_seq=seq[Lbar_pos:(Lbar_pos+len(Lbar))]
	Bar_seq=seq[Bar_pos:(Bar_pos+Nbar)]
	Rbar_seq=seq[Rbar_pos:(Rbar_pos+len(Rbar))]
	if lv.distance(Lbar,Lbar_seq)>Lbar_dist:
		return 'Lbar_incorrect',0 #Lbar sequence was too different from expected
	if lv.distance(Rbar,Rbar_seq)>Rbar_dist:
		return 'Rbar_incorrect',0 #Rbar sequence was too different from expected
	if seq[(Bar_pos-Bar_clamp):Bar_pos]!=Lbar[-Bar_clamp:]:
		return 'Lclamp_mismatch',0 #clamp left of barcode mismatch
	if seq[Rbar_pos:(Rbar_pos+Bar_clamp)]!=Rbar[0:Bar_clamp]:
		return 'Rclamp_mismatch',0 #clamp right of barcode mismatch
#	return Bar_seq,(Rbar_pos+len(Rbar))
	return Bar_seq,Rbar_pos #don't trim Rbar


def Radapter_trim_fuzzy(seq,Lbar,Nbar=20):
#takes in sequence in string form, outputs index at which trimming should occur and barcode sequence
    matches_l=fs.find_near_matches(Lbar,seq,max_l_dist=4)
    if matches_l==[]:
        return 'Lbar_not_found',-1
    else:
        dist_l=[i.dist for i in matches_l]
        min_index_l=dist_l.index(min(dist_l))
        match_l=matches_l[min_index_l]
    
    #use the min function in case we have a very short oligo and sequence pass the barcode
    #use code below for ignoring barcodes not Nbar length
    #but we can keep them, since we only use barcodes to link oligos to barcodes
    #even if a barcode is 18, 19bp, it may still be used
    #N=min([len(seq)-match_l.end, Nbar])
    #if N!=Nbar:
    #	return 'Barcode_incorrect_length',-1
    #else:
    #	bar=seq[match_l.end:match_l.end+Nbar]
    #	trim_index=match_l.start
    #again, note in this case, the barcode may not be exactly 20bp long
    bar=seq[match_l.end:match_l.end+Nbar]
    trim_index=match_l.start		
    return bar,int(trim_index)


def Ladapter_trim_fuzzy(seq,Ladp):
#takes in sequence in string form, outputs index at which trimming should occur
	matches=fs.find_near_matches(Ladp,seq,max_l_dist=3)
	if matches==[]:
		return len(seq)
	else:
		dist=[i.dist for i in matches]
		min_index=dist.index(min(dist))
		match=matches[min_index]
		trim_index=match.start
		
	return int(trim_index)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_file_path',type=str)
    parser.add_argument('-o',type=str,help='output file path for trimmed sequences')
    parser.add_argument('-d',type=str,help='output file path for discarded sequences')
    parser.add_argument('-Lbar',default='GATCGCGTCGACGAACCTCTAGA',type=str,help='sequence left of barcode') #MPRADel
    #it is the read 1 read primer, no need to include it
    #parser.add_argument('-Rbar',default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',type=str,help='sequence right of barcode') 
    parser.add_argument('-Ladp',default='ACTGGCCGCTTGACG',type=str,help='sequence of right adapter')
    parser.add_argument('-barxN',action='store_true')
    parser.add_argument('-seqxN',action='store_true')
    #whether or not to trim R adapter by just bp 
    parser.add_argument('-Ladp_trim_bp',action='store_true')
    
    parsed = parser.parse_args()
    #hard code for now
    Ladp_trim_bp_cutoff=16

    fastq_file = gzip.open(parsed.fastq_file_path,'rb')
    trimmed_file = gzip.open(parsed.o,'wb')
    discard_file = gzip.open(parsed.d,'wb')
 
    #	it=0
    #	Nline=0
    for record in SeqIO.parse(fastq_file,'fastq'):
        
        #rev complement sequence first
        description=record.description
        #reverse complement sequence, not it is important to store description beforehand as that info is lost after reverse complement
        #https://www.biostars.org/p/303826/
        record=record.reverse_complement()
        seq=str(record.seq)
        
        if parsed.seqxN: #if the seqxN option is used, skip sequences with N's anywhere.
            if 'N' in seq:
                description = description + ' N_in_seq' #add barcode to description
                record.description='_'.join(description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
                record.id=record.description #so separate fields aren't added
                record.name=record.description #so separate fields aren't added
                SeqIO.write([record],discard_file,'fastq')
                continue

        [barcode,r_index] = Radapter_trim_fuzzy(seq,parsed.Lbar)
        if parsed.barxN: #if the barxN option is used, skip sequences with N in the barcode position.
            if 'N' in barcode:
                description = description + ' bar_N:%s' %(barcode) #add barcode to description
                record.description='_'.join(description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
                record.id=record.description #so separate fields aren't added
                record.name=record.description #so separate fields aren't added
                SeqIO.write([record],discard_file,'fastq')
                continue

 
        if r_index==-1: #failed to identify the barcode sequence
            description = description + ' ' + barcode
            record.description='_'.join(description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
            record.id=record.description #so separate fields aren't added
            record.name=record.description #so separate fields aren't added 
            SeqIO.write([record],discard_file,'fastq')
            continue
       
        #good sequence 
        description = description + ' bar:%s' %(barcode) #add barcode to description
        record.description='_'.join(description.split(' ')) #convert ' ' to '_' so samtools doesn't cut them off
        record.id=record.description #so separate fields aren't added
        record.name=record.description #so separate fields aren't added
      
        if parsed.Ladp_trim_bp:
            #simply trim off first Ladp_trim_bp_cutoff bp, we do this instead of fuzzysearch due to first several bp being especially noisy
            l_index=Ladp_trim_bp_cutoff
        else:
            l_index = Ladapter_trim_fuzzy(seq,parsed.Ladp)
        
        record=record[l_index:r_index]
    #		record=record[r_index:] #don't trim Ladapter

        SeqIO.write([record],trimmed_file,'fastq')
    #		trimmed.append(record)

    fastq_file.close()
    trimmed_file.close()
    discard_file.close()


def one_record():
#	Lbar='CGAGCTGTACAAGTAATTCTAGTTG' #PGKGiga
	Lbar='GTTTAAAGCCCAACGCTAGTC' #BmtGiga
	Rbar='CGAGCTCGCTAGCCT'
	Ladp='AGATCGGAAGAGCGTCG'
	
	fastq_file=open('/cluster_path/sequencing/AALE9_2/fastq/AALE9.AACTTGAC.extendedFrags.fastq','r')
	record=SeqIO.parse(fastq_file,'fastq').next()

	seq=str(record.seq)
		
	[barcode,r_index] = Radapter_trim(seq,Lbar,Rbar)
	if r_index==0: #failed to identify the barcode sequence
		print 'DISCARDED SEQUENCE'
		print record.format('fastq')
		fastq_file.close()
		return
	record.description = record.description + ' bar:%s' %(barcode)
	record.description='_'.join(record.description.split(' '))
	record.id=record.description
	record.name=record.description
	
	l_index = Ladapter_trim(seq,Ladp)
	record=record[l_index:r_index]
#	record=record[r_index:] #don't trim Ladapter

	print 'TRIMMED SEQUENCE'
	print record.format('fastq')
	
	fastq_file.close()


if __name__ == '__main__':
	if len(sys.argv)==1:
		one_record()
	else:
		main()
