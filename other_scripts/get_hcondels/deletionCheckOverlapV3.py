from __future__ import with_statement
import numpy as np
import sys

#this script takes in an intersected file of a region with extended coordinates with the VCF annotated deletion
#then checks to see if the deletions are removing the same site.
#to see this, we simply remove the specified deletion sequence and check to see if 
#the resulting sequence is the same

arg=sys.argv
bedIntersectedFileName=arg[1]
deletionSequenceFileName=arg[2]

#hybOverStartInd is where the start index is for the region over our site of interest
#hybUnderStartInd is where the start index is for the region under our site of interest
#VCFStartInd is the position of the VCF start position

hybOverStartInd=int(arg[3])
hybUnderStartInd=int(arg[4])
VCFStartInd=int(arg[5])

#invInd is whether or not the sequence is inverted
invInd=VCFStartInd-2
#this is a simple file of chimp id in first column and +/- in second column indicating whether or not chimp sequence is inverted in human sequence
outFileName=arg[6]

#read in chimp fasta sequence
chimp_fasta_seq=dict()
with open(deletionSequenceFileName, "r") as deletionSequenceFile:
    for ind, l1 in enumerate(deletionSequenceFile):
        line=l1.split()
        lineEntry=line[0]
        #should be alternating between loading chr into dictionary and putting sequence into dictionary
        #each entry is a coordinate 
        if(">" == lineEntry[0]):
            chr_name=line[0][1:len(lineEntry)]
            chimp_fasta_seq[chr_name]=[]
        else:
            chimp_fasta_seq[chr_name]=list(lineEntry)

allLines=[]
with open(bedIntersectedFileName, "r") as bedIntersectedFile:
    for l1 in bedIntersectedFile:
        line=l1.split()
        
        chrName=line[0]
            
        chimpHumanHybridInterestOverStartCoord=int(line[hybOverStartInd])
        chimpHumanHybridInterestOverEndCoord=int(line[hybOverStartInd+1])
   
        coordId='|'.join([chrName, str(chimpHumanHybridInterestOverStartCoord), str(chimpHumanHybridInterestOverEndCoord)])
 
        chimpHumanHybridInterestUnderStartCoord=int(line[hybUnderStartInd])
        chimpHumanHybridInterestUnderEndCoord=int(line[hybUnderStartInd+1])

        chimpHumanHybridInterestUnderVCFStartCoord=int(line[VCFStartInd])
        chimpHumanHybridInterestUnderVCFEndCoord=int(line[VCFStartInd+1])

        startPosVCF=chimpHumanHybridInterestUnderVCFStartCoord-chimpHumanHybridInterestOverStartCoord
        endPosVCF=chimpHumanHybridInterestUnderVCFEndCoord-chimpHumanHybridInterestOverStartCoord
        
        startPosOrig=chimpHumanHybridInterestUnderStartCoord-chimpHumanHybridInterestOverStartCoord
        endPosOrig=chimpHumanHybridInterestUnderEndCoord-chimpHumanHybridInterestOverStartCoord

        invInfo=line[invInd]
        chimp_seq_hybrid=chimp_fasta_seq[coordId]
        
        #remove deleted sequence by slicing
        #remember that we are in BED format, so include first position, but don't include last position
        #remember VCF coordinate includes first bp which is not part of deletion, so need to add 1 to include that base
        
        chimp_vcf_seq_deletion_removed=chimp_seq_hybrid[0:startPosVCF+1]+chimp_seq_hybrid[endPosVCF:len(chimp_seq_hybrid)]
        chimp_orig_seq_deletion_removed=chimp_seq_hybrid[0:startPosOrig]+chimp_seq_hybrid[endPosOrig:len(chimp_seq_hybrid)]
        #print(chimp_seq)
        #print(chimp_vcf_seq_deletion_removed)
        #print(chimp_orig_seq_deletion_removed) 
        
        #compare vcf deleted sequence with original deleted sequence
        seqMatch=0
        #reverse sequence for proper comparison
        if(chimp_vcf_seq_deletion_removed==chimp_orig_seq_deletion_removed):
            #for debugging
            #print(line)
            #print(''.join(chimp_seq_hybrid))
            #print(''.join(chimp_seq_hybrid[startPosOrig:endPosOrig]))
            #print(''.join(chimp_seq_hybrid[startPosVCF+1:endPosVCF]))
            seqMatch=1            
        
        newLine=line+[seqMatch]
        allLines.append(newLine)
      
#last field in file will be whether or not there is a sequence match or not
#simply write out whether or not VCF file annotated deletion site matches 
#can simply write out original file with an extra field of matching or not
with open(outFileName, "w") as outFile:
    for ind in range(len(allLines)):
        for ind2 in range(len(allLines[ind])):
            outFile.write(str(allLines[ind][ind2])+"\t")   
        outFile.write("\n")



    



