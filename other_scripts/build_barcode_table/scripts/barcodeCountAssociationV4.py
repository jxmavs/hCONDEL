import sys
import pysam
import gzip

arg=sys.argv

#link the read to the barcode
#read mapped to library file list
bamListFileName=arg[1]
barcodeInfoOutFileName=arg[2]
seqNamesFileName=arg[3]

#bamListFileName="/cluster_path_temp/nova_12-15-18/mpradel-ehb/enhancer_barcode_1_22_19_merged_combined_file_list_bam.aligned.txt"

def findFlag(lineArr, flagStr):
    for ind in range(len(lineArr)):
        if(flagStr == lineArr[ind][0]):
            return ind


def getMisAlignmentStats(alignInfo):    
    XMind=findFlag(alignInfo, "XM")
    XMInfo=alignInfo[XMind][1]
    numMismatch=int(XMInfo)
    
    XOind=findFlag(alignInfo, "XO")
    XOInfo=alignInfo[XOind][1]
    numGaps=int(XOInfo)
    #print(XOFlag)
    
    XGind=findFlag(alignInfo, "XG")
    XGInfo=alignInfo[XGind][1]
    numGapExtensions=int(XGInfo)
        
    allMisStats=[numMismatch, numGaps, numGapExtensions]    
    return(allMisStats)


#store seqNamesDict to allow less memory usage in following code
#seqNamesArr is for mapping ids back to sequences for writing out
seqNamesDict=dict()
seqNamesArr=[]
i=0
with open(seqNamesFileName) as seqNamesFile:
    for line in seqNamesFile:
        line=line.split()
        seqNamesDict[line[0]]=i
        seqNamesArr.append(line[0])
        i+=1

#keep track of barcodes that map to multiple sequences
barcodeArr=[]
barcodeSeqDict=dict()
barcodeArrLength=0
print("begin")
#read in list of bam files
with open(bamListFileName, "r") as bamListFile:
    
    for line0 in bamListFile:
        
        bamFileName=line0.split()[0]
        bamFile=pysam.AlignmentFile(bamFileName, "rb")

        print(bamFileName)
        print(barcodeArrLength)
        #flush out output
        sys.stdout.flush()         
        for line in bamFile.fetch(until_eof=True):
             
            readName=line.query_name
            seqName=line.reference_name
            seqId=seqNamesDict[seqName]
            #get barcode, split instead to get barcodes that are less than 20 bp
            readNameParts=readName.split(":")
            barcode=readNameParts[len(readNameParts)-1]

            readLen=line.query_length
            #first entry is number of reads
            allStats=[1, readLen]
            allStats.extend(getMisAlignmentStats(line.tags))
            
            #new entry
            if barcode not in barcodeSeqDict:
                #dupInfo stores whether or not information is a duplicate            
                dupInfo=0
                #the second entry stores index of where barcodeSeqDict is
                barcodeSeqDict[barcode]=[seqId, barcodeArrLength]
                barcodeArr.append([allStats, [seqId], dupInfo])
                barcodeArrLength+=1
                #sys.stdout.flush()
                #print(barcodeArrLength)
            else:
                index=barcodeSeqDict[barcode][1]
               
                #initial check for efficiency purposes 
                if barcodeSeqDict[barcode][0]!=seqId:
                    #add on sequence, it is a duplicate
                    if seqId not in barcodeArr[index][1]:
                        barcodeArr[index][1].append(seqId)
                        #add on number of duplicates barcode has seen
                        barcodeArr[index][2]+=1
                    
                #sum up total stats
                for ind in range(len(barcodeArr[index][0])):
                    barcodeArr[index][0][ind]+=allStats[ind]
        
        print(barcodeArrLength) 
        sys.stdout.flush()
        bamFile.close()

        
#the output is:
#barcode, number of total sequences linked to barcode, total number of bases linked to barcode, total number of mismatches linked to barcode, total number of gaps linked to barcode, total size of gaps in bp linked to barcode, ids of sequences associated to barcode seperated by ":", whether or not barcode is linked to multiple DIFFERENT sequences 
with open(barcodeInfoOutFileName, "w") as barcodeOutFile:
    for barcode in barcodeSeqDict.keys():
        index=barcodeSeqDict[barcode][1]
        info=barcodeArr[index]
        barcodeOutFile.write(barcode+"\t")
        
        #write out stats
        for ind in range(len(info[0])):
            barcodeOutFile.write(str(info[0][ind])+"\t")

        seqStr=""
        for ind in range(len(info[1])-1):
            seqStr+=str(seqNamesArr[info[1][ind]])+":"
        seqStr+=str(seqNamesArr[info[1][len(info[1])-1]])
        
        #write out the sequence associatd with barcode
        barcodeOutFile.write(seqStr+"\t")

        #write out number of duplicates associated with barcode
        barcodeOutFile.write(str(info[2]))
        barcodeOutFile.write("\n")

            
