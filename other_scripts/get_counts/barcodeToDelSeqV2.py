import sys
import gzip
import fuzzysearch as fs
import Levenshtein as lv
from Bio import SeqIO
from Bio.Seq import Seq

#this script takes a fastq file and links the barcode to deletion sequence by previous enhance/barcode table
# a table of counts is produced per del sequence
arg=sys.argv

#link the read to the barcode
barcodeFastqFileName=arg[1]
barcodeAssocTableFileName=arg[2]

#seqName file is just a file of the names of all deletion sequences
seqNamesFileName=arg[3]
barcodeLength=int(arg[4])
seqIdFilterPct=float(arg[5])

outFileName=arg[6]


#I can use ATCGCGTCGACGAACCTCTAGA
#I can perhaps use GACGAACCTCTAGA instead and max_l_dist of 2
def Radapter_trim_fuzzy(seq,Lbar="AATCGCGTCGACGAACCTCTAGA",Nbar=20):
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
    #   return 'Barcode_incorrect_length',-1
    #else:
    #   bar=seq[match_l.end:match_l.end+Nbar]
    #   trim_index=match_l.start
    #again, note in this case, the barcode may not be exactly 20bp long
    bar=seq[match_l.end:match_l.end+Nbar]
    trim_index=match_l.start
    return bar,int(trim_index)



#set the count of all initial sequences to 0
seqCountDict=dict()
with open(seqNamesFileName, "r") as seqNamesFile:
    for line in seqNamesFile:
        seqName=line.split()[0]
        seqCountDict[seqName]=[]

badDupBarcodeDict=dict()
badSynthErrorBarcodeDict=dict()
barcodeDict=dict()
with open(barcodeAssocTableFileName, "r") as barcodeAssocTableFile:
    for line in barcodeAssocTableFile:
        line=line.split('\t')
        barcode=line[0]
        seqName=line[6]
        dupInfo=int(line[len(line)-1])
        #skip duplicated sequences 
        if dupInfo>0:
            badDupBarcodeDict[barcode]=0
            continue

        #skip bad barcodes 
        barcodeArrStats=[ int(line[ind]) for ind in range(1, 6) ]

        totalReadLen=barcodeArrStats[1]
        numMismatches=barcodeArrStats[2]
        numGaps=barcodeArrStats[4]
        
        if float(numGaps+numMismatches)/totalReadLen > 1-seqIdFilterPct:
            badSynthErrorBarcodeDict[barcode]=0
            continue
        
        barcodeDict[barcode]=seqName
        
numBarcodesNotLinked=0
totalNumberofBarcodes=0
totalNumberSeq=0
numSeqNoAdapter=0
numBarcodesLinkedToBadDupBarcode=0
numBarcodesLinkedToBadSynthErrorBarcode=0

with gzip.open(barcodeFastqFileName, "rb") as barcodeFastqFile:
    
    for record in SeqIO.parse(barcodeFastqFile,'fastq'):
        record=record.reverse_complement()
        seq=str(record.seq)        

        [barcode,trim_index]=Radapter_trim_fuzzy(seq)
        
        #barcode=line[0][0:barcodeLength]
        totalNumberSeq+=1
        if barcode=='Lbar_not_found':
            numSeqNoAdapter+=1
        elif barcode in badDupBarcodeDict:
            totalNumberofBarcodes+=1
            numBarcodesLinkedToBadDupBarcode+=1
        elif barcode in badSynthErrorBarcodeDict:
            totalNumberofBarcodes+=1
            numBarcodesLinkedToBadSynthErrorBarcode+=1
        elif barcode not in barcodeDict:
            totalNumberofBarcodes+=1
            numBarcodesNotLinked+=1
        else:
            totalNumberofBarcodes+=1
            seqName=barcodeDict[barcode]
            seqCountDict[seqName].append(barcode)

with open(outFileName+".stats", "w") as outFileStats:
    outFileStats.write("Num_Total_Seq\tNum_Total_Barcodes\tNum_Barcodes_Not_Linked\tNum_Bad_Dup_Barcodes_Linked\tNum_Bad_Synth_Error_Barcodes_Linked\tNum_Seq_No_Adapter\n")
    outFileStats.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(totalNumberSeq, totalNumberofBarcodes, numBarcodesNotLinked, numBarcodesLinkedToBadDupBarcode, numBarcodesLinkedToBadSynthErrorBarcode, numSeqNoAdapter))

#numBarcodesLinkedToBadDupBarcode+numBarcodesLinkedToBadSynthErrorBarcode is the total number of barcodes linked to bad barcodes
#print("Number of total barcodes: {0}\nNumber of barcodes that did not match: {1}\n".format(totalNumberofBarcodes, numBarcodesNotLinked))
#out file is three columns - one is column name, the second is the count of unique barcodes, the third is the total count
with open(outFileName, "w") as outFile:
    for seq in sorted(seqCountDict.keys()):
        info=seqCountDict[seq]
        uniqBarcodeCount=len(set(info))
        totalCount=len(info)
        outFile.write(seq+"\t"+str(uniqBarcodeCount)+"\t"+str(totalCount)+"\n")

        
