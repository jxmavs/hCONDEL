from __future__ import with_statement
import numpy as np
import sys

arg = sys.argv

#input contains two coordinates, the conserved coordinate to obtain the conserved sequence and the deletion coordinate to check that the site is deleted everywhere
inputFileName=arg[1]
#mafdir contains the directory with each file containg a MAF sequence
mafDir=arg[2]
outFileName=arg[3]
#check if there is sequence present, as well as if the sequence matches exactly to other sequences
#output should be one info for all primates, as well as one output for rest of mammals
#I need to check the flanking MSA region as well, not just the single bp

#the categories below are for checking the underlying bases/sequences aligning to each specific category
chimp_macaque_ids=["panTro4", "rheMac8"]
chimp_macaque_mouse_ids=["panTro4", "rheMac8", "mm10"]
all_primate_ids=["panTro4", "panPan1", "gorGor4", "ponAbe2", "rheMac8"]
all_species_ids=["panTro4", "panPan1", "gorGor4", "ponAbe2", "rheMac8", "mm10", "bosTau8", "canFam3", "monDom5", "galGal4", "ornAna2"]
all_categories=[chimp_macaque_ids,chimp_macaque_mouse_ids, all_primate_ids, all_species_ids]

        

def getOverlapInd(list1, list2):
    allIndKeep=[]
    for ind in range(len(list1)):
        if(list1[ind] in list2):
            allIndKeep.append(ind)
    return(allIndKeep)

#this function calculates the proportion of A,T,G,C as well as G/C (G+C) 
def getSeqContent(seq_arr):
    allSeqLen=0
    allSeqContent=[0]*4
    #seqInterest=[["A","a"],["T", "t"], ["C", "c"], ["G","g"]]
    seqIndeces=dict()
    seqIndeces["A"]=0
    seqIndeces["T"]=1
    seqIndeces["G"]=2
    seqIndeces["C"]=3
    seqIndeces["a"]=0
    seqIndeces["t"]=1
    seqIndeces["g"]=2
    seqIndeces["c"]=3
    for ind in range(len(seq_arr)):
        #ignore gaps
        #allSeqLen+=sum(seq_arr[ind,:]!="-")
        for ind2 in range(len(seq_arr[ind])):
            #iterate over base
            for ind3 in range(len(seq_arr[ind][ind2])):
                if(seq_arr[ind][ind2][ind3] in ["A","a", "T", "t", "C", "c", "G","g"]):
                    
                    indChar=seqIndeces[seq_arr[ind][ind2][ind3]]
                    allSeqContent[indChar]+=1
                    allSeqLen+=1
                
            #seqInterestArr=seqInterest[ind2]
            #allSeqContent[ind2]+=sum([i in seqInterestArr for i in seq_arr[ind,:] ])

    #calculate percentages    
    seqPctArr=[0]*5
    for ind in range(len(allSeqContent)):
        seqPctArr[ind]=allSeqContent[ind]/float(allSeqLen)
    #below is to calculate overall GC content
    seqPctArr[len(seqPctArr)-1]=seqPctArr[len(seqPctArr)-2]+seqPctArr[len(seqPctArr)-3]
    return(seqPctArr)

#this function reads in a MAF file and concatenates MAF sequences together
#it returns the MAF exact sequences, the species in the sequence,  the start/end positions, and finally whether or not the MAF has gaps between blocks (I have not seen this yet)
def readMAFFile(mafFile):
    with open(mafFile, "r") as idFile:
        #ignore first line
        idFile.readline()
        allMAFSeq=[]
        allMAFSpecies=[]
        blockStartEndPosArr=[]

        MAFSeqIter=[]
        MAFSeqSpeciesIter=[]
        blockStartEndPosArrIter=[]

        #blockInd keeps track of which MAF block we're looking at
        blockInd=-1

        for l1 in idFile:
            line=l1.split()
            #print(line)
            #if we have an empty line, then its a blank line
            if(len(line)==0):
                continue
            elif("score" in line[1]):
                blockInd+=1
                #print("score") 
                if(blockInd>0):
    
                    allMAFSeq.append(MAFSeqIter)
                    allMAFSpecies.append(MAFSeqSpeciesIter)
                    blockStartEndPosArr.append(blockStartEndPosArrIter)
                
                    MAFSeqIter=[]
                    MAFSeqSpeciesIter=[]
                    blockStartEndPosArrIter=[] 
            else:

                seqAlignment=line[6]
                
                if(seqAlignment=="(null)"):
                    continue

                seqAlignChar=list(seqAlignment)
                species=line[1].split(".")[0]
                #print(species) 
                #bed format positions for each block, keep track of them
                seqAlignStartPos=int(line[2])
                seqAlignEndPos=seqAlignStartPos+int(line[3])
                blockStartEndPosArrIter.append([seqAlignStartPos, seqAlignEndPos])

                
                #seqAlignment.replace("-","")  
                MAFSeqIter.append(seqAlignChar)
                MAFSeqSpeciesIter.append(species)

        #add on last block

        allMAFSeq.append(MAFSeqIter)
        allMAFSpecies.append(MAFSeqSpeciesIter)
        blockStartEndPosArr.append(blockStartEndPosArrIter)
        

        #print("tot blocks")
        #print(blockInd)
        
        #initially consider all species in window as "common"
        
        #get common species for all MAFs - for all MAF blocks, we take the intersection of the common species
        #initialize by setting all common species to the species in the first block
        common_species_arr=allMAFSpecies[0]
    
        for ind in range(1, len(allMAFSpecies)):
            #iterate through all species 
            common_species_ind_arr=[0]*(len(common_species_arr))
            for ind1 in range(len(allMAFSpecies[ind])):
                species=allMAFSpecies[ind][ind1]

                foundInd=-1
                for i, j in enumerate(common_species_arr):
                    if(j==species):
                        foundInd=i
                        break

                if(foundInd!=-1):
                    common_species_ind_arr[foundInd]=1
            
            #keep species that do overlap, remove species that don't overlap
            common_species_temp_arr=[]
            for ind1 in range(len(common_species_ind_arr)):
                if(common_species_ind_arr[ind1]==1):
                    common_species_temp_arr.append(common_species_arr[ind1])
            common_species_arr=common_species_temp_arr
       
        #iterate through the categories 
        #read in initial sequence       
        seq_iter_arr=[] 
        for ind1 in range(len(common_species_arr)):
            commonSpecies=common_species_arr[ind1]
            indKeep=allMAFSpecies[0].index(commonSpecies)
            #print("indkeep")
            #print(indKeep)
            #print(allMAFSeq[0])
            seq_iter_arr.append(allMAFSeq[0][indKeep])
        all_seq_arr=np.array(seq_iter_arr)
        #this should always return a 0?
        chimpInd=allMAFSpecies[0].index("panTro4")
        
        prevBlockEndPos=blockStartEndPosArr[0][chimpInd][1]
        #print(beginBlockInd)
        #print(endBlockInd)
        numMAFGaps=0

        #print(common_species_arr)
        for ind in range(1, len(blockStartEndPosArr)):
            seq_iter_arr=[]
            #search through to find common sequences, could be possible that its not in the same order?
            #print(allMAFSeq[beginBlockInd])
            for ind1 in range(len(common_species_arr)):
                commonSpecies=common_species_arr[ind1]
                indKeep=allMAFSpecies[ind].index(commonSpecies)
                #print(indKeep)
                seq_iter_arr.append(allMAFSeq[ind][indKeep])

            #merge into common MAF sequence, combine via columns 
            #calculate difference between position between previous MAF block and current MAF block
            #if there is a gap between MAF's (I don't expect there to be any)  make a column for the gaps and keep track of it
            chimpInd=allMAFSpecies[ind].index("panTro4")
            currentBlockStartPos=blockStartEndPosArr[ind][chimpInd][0]
            #print(currentBlockStartPos)
            if(currentBlockStartPos-prevBlockEndPos>0):
                gapLen=currentBlockStartPos-prevBlockEndPos
                numMAFGaps+=gapLen
                gapArr=np.full((len(common_species_arr), gapLen), "-", dtype="string")
                all_seq_arr=np.hstack((all_seq_arr, gapArr))
                all_seq_arr=np.hstack((all_seq_arr, seq_iter_arr))
                print("MAFGap!")
            else:
                all_seq_arr=np.hstack((all_seq_arr, seq_iter_arr))            
            #merge into common MAF sequence, combine via columns 
            #calculate difference between position between previous MAF block and current MAF block
            #if there is a gap between MAF's (I don't expect there to be any)  make a column for the gaps and keep track of it
            chimpInd=allMAFSpecies[ind].index("panTro4")
            currentBlockStartPos=blockStartEndPosArr[ind][chimpInd][0]
            #print(currentBlockStartPos)
            all_seq_arr=np.hstack((all_seq_arr, seq_iter_arr))

            prevBlockEndPos=blockStartEndPosArr[ind][chimpInd][1]
    #blockStartEndPosArr[beginBlockInd:(endBlockInd+1)] contains the block positions
    return([ all_seq_arr, common_species_arr, blockStartEndPosArr, numMAFGaps ])
 
coordArr=[]
allPosInfo=[]
allPosMetaInfo=[]
with open(inputFileName , "r") as inputFile:
    #enumerate through input file , each line is a coordinate
    for ind, l1 in enumerate(inputFile):
        l1=l1.split()
        consCoordChr=l1[0]
        consCoordStart=int(l1[1])
        consCoordEnd=int(l1[2])
        
        coordArr.append(l1)
        idFileName=mafDir+"/"+"|".join(l1[0:3])+".maf"

        #print(idFileName)
        #read in MAF file first and obtain sequences, species, block start and end positions
        all_seq_info=readMAFFile(idFileName)
        all_seq_np=all_seq_info[0]
        all_seq_species=all_seq_info[1]
        all_seq_block_coord=all_seq_info[2]
        numMAFGaps=all_seq_info[3] 
        allCategoryInfo=[]


        for ind in range(len(all_categories)):
            category=all_categories[ind]
            #pull out only sequences in the specific category
            overlapInd=getOverlapInd(all_seq_species, category)
            
            #get the number of species in each category
            numSpeciesInCategory=len(overlapInd)
        
            #get the actual species in the form of an ID
            allSpeciesInCat="|".join([ all_seq_species[j] for i,j in enumerate(overlapInd) ])
       
            #get A,T,G,C, G/C content for sequencees in the category
            seqContentArr=getSeqContent(all_seq_np[overlapInd,]) 
            seqContentArr.extend([numSpeciesInCategory, allSpeciesInCat])
            allCategoryInfo.append(seqContentArr)
 
        #metaInfo is independent of categories
        allPosMetaInfo.append([numMAFGaps])
        allPosInfo.append(allCategoryInfo)


with open(outFileName, "w") as outFile:
    for ind in range(len(allPosInfo)):
        #write out the original coordinates
        for ind1 in range(len(coordArr[ind])):
            outFile.write(str(coordArr[ind][ind1])+"\t")

        #number of iterations ofthe first loop in the loop below is the number of categories
        for ind2 in range(len(allPosInfo[ind])):
            for ind3 in range(len(allPosInfo[ind][ind2])):
                outFile.write(str(allPosInfo[ind][ind2][ind3])+"\t")

        #write out the meta info, includes the length of the deletion, actual window length used (since the window length could be truncated), the number of gaps between MAF blocks for the window position
        for ind3 in range(len(allPosMetaInfo[ind])):
            outFile.write(str(allPosMetaInfo[ind][ind3])+"\t")
        outFile.write("\n")



