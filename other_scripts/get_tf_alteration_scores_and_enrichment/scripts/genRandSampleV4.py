import sys
import numpy as np
#this script takes a raw sam aligned file and checks to see the number of mismatches/gaps for readName==refName
#genRandSampleV2.py takes in file with human/chimp differences and draws based on that

arg = sys.argv

origSuperSetFileName=arg[1]
origFileName=arg[2]
nrowOrigFile=int(arg[3])
lengthBoundPct=float(arg[4])
totalMisMatchBoundPct=float(arg[5])
GCBoundPct=float(arg[6])
logOddsBoundPct=float(arg[7])
sampledOutFileName=arg[8]

#origSuperSetFileName="/cluster_path_temp/roadmap/LCL_data/mostConservedMapped_new_0.001_human_size_filtered_extended_200_hg38_aligned_panTro4_check_stats_revised.txt"
#origFileName="/cluster_path_temp/roadmap/LCL_data/orig_cons_af_1_final_filtered_extended_200_hg38_aligned_panTro4_check_stats_revised.txt"
#nrowOrigFile=11838
#lengthBoundPct=0.01
#totalMisMatchBoundPct=0.01
#sampledOutFileName="/cluster_path_temp/roadmap/LCL_data/test_test_cons_sequences_info.txt"
#with open(origSuperSetFileName, "r") as origSuperSetFile:
#    for ind, l1 in enumerate(origSuperSetFile):
#        line=l1.split()
#        seqID=line[0]
#        numTotalMisMatch=line[1]

superSetData=np.genfromtxt(origSuperSetFileName, dtype=None, delimiter="\t")
superSetData.dtype.names=('seqId', 'chrID','consSeqLen', 'totalMisMatchPct', 'GCPct', 'logOdds')

allRandSeq=np.empty([nrowOrigFile ,6], dtype=object)

#seqLowBound=4.95
#seqUpBound=5.05
#totalMisMatchLowBound=0
#totalMisMatchUpBound=100
#read in file of conserved blocks with deletions
#keep track of all seqIDs read 
allSampledSeqIDs=dict()
with open(origFileName, "r") as origFile:
    for ind, l1 in enumerate(origFile):
        #print(ind)
        line=l1.split()
        chrID=line[1]
        seqLen=int(line[2])
        totalMisMatchPct=float(line[3])
        GCPct=float(line[4])
        logOdds=int(line[5])
        #the first 6 coordinates are the original conserved sequences
        consID='|'.join(line[0].split("|")[0:6])

        #if the conserved sequence is a duplicate, then this means that we have multiple deletions lieing in the sequence
        if(consID in allSampledSeqIDs):
            indOrig=allSampledSeqIDs[consID]
            allRandSeq[ind,:]=allRandSeq[indOrig,:]    
        else:
            subsetInd=np.array([])

            #reset to original threshold values
            lengthBoundPctIter=lengthBoundPct
            totalMisMatchBoundPctIter=totalMisMatchBoundPct
            GCBoundPctIter=GCBoundPct
            logOddsBoundPctIter=logOddsBoundPct
            #if the length of the array is zero, then reget the indeces by relaxing the cutoffs
            #otherwise keep the subset,
            while True:
         
                seqLowBound=max([0, seqLen-lengthBoundPctIter*(seqLen)])
                seqUpBound=seqLen+lengthBoundPctIter*(seqLen)

                #totalMisMatchLowBound=max([0, totalMisMatchPct-totalMisMatchBoundPctIter*(totalMisMatchPct)])
                #totalMisMatchUpBound=totalMisMatchPct+totalMisMatchBoundPctIter*(totalMisMatchPct)
                
                totalMisMatchLowBound=max([0, totalMisMatchPct-totalMisMatchBoundPctIter])
                totalMisMatchUpBound=totalMisMatchPct+totalMisMatchBoundPctIter
                
                #GCLowBound=max([0, GCPct-GCBoundPctIter*(GCPct)])
                #GCUpBound=GPct+GCBoundPctIter*(GCPct)

                GCLowBound=max([0, GCPct-GCBoundPctIter])
                GCUpBound=GCPct+GCBoundPctIter
               
                logOddsLowBound=max([0, logOdds-logOddsBoundPctIter*(logOdds)])
                logOddsUpBound=logOdds+logOddsBoundPctIter*(logOdds)
                
                #print(seqLowBound)
                #print(seqUpBound)
                #print(totalMisMatchLowBound)
                #print(totalMisMatchUpBound)
                
                subsetInd=( (superSetData['chrID']==chrID) & (superSetData['consSeqLen']>=seqLowBound) & (superSetData['consSeqLen']<=seqUpBound) & (superSetData['totalMisMatchPct']>=totalMisMatchLowBound) & (superSetData['totalMisMatchPct']<=totalMisMatchUpBound) & (superSetData['GCPct']>=GCLowBound) & (superSetData['GCPct']<=GCUpBound) & (superSetData['logOdds']>=logOddsLowBound) & (superSetData['logOdds']<=logOddsUpBound) ).nonzero()[0]
                
                #subsetInd=( (superSetData['consSeqLen']>=seqLowBound) & (superSetData['consSeqLen']<=seqUpBound) ).nonzero()[0] 
                #check below is to see if there are no elements in the subset
                #theoretically, it should never be a 0 due to the conserved sequence file being a superset of the orignal file
                #but can be 0 if there is differences in alignment stats between runs of bowtie
                #if(len(subsetInd)==0):
                #    print(line)
                

                if(len(subsetInd)>1):
                    break
                else:
                    #increment by 0.01 if we fail to randomly sample anything
                    totalMisMatchBoundPctIter=totalMisMatchBoundPctIter+0.01
                    lengthBoundPctIter=lengthBoundPctIter+0.05
                    GCBoundPctIter=GCBoundPctIter+0.03
                    logOddsBoundPctIter=logOddsBoundPctIter+0.05
                    #print(str(len(subsetInd))+"\t"+str(GCBoundPctIter))
                    #print(line)
                #randomly sample sequnce which matches profile


            randInd=np.random.randint(0, len(subsetInd))     
            sampleInd=subsetInd[randInd]
 
            allRandSeq[ind,:]=[ superSetData['seqId'][sampleInd], len(subsetInd), lengthBoundPctIter, totalMisMatchBoundPctIter, GCBoundPctIter, logOddsBoundPctIter] 
            #print(allRandSeq[ind,:])
            #after sampling, remove sequence from being sampled again
            superSetData=np.delete(superSetData, sampleInd, 0) 

            allSampledSeqIDs[consID]=ind

#np.savetxt(sampledOutFileName, allRandSeq, delimiter='\t')
with open(sampledOutFileName, "w") as sampledOutFile:
    for ind in range(len(allRandSeq)):
        for ind2 in range(len(allRandSeq[ind])):
            sampledOutFile.write(str(allRandSeq[ind][ind2])+"\t")
        sampledOutFile.write("\n") 
#for each conserved block that has a deletion, sample randomly from the superset file
#get all blocks that match the profile, then from those blocks, sample randomly






