import sys
#this script takes the table of barcode counts and converts it to a count per del
#it removes barcodes which are mapped to multiple del sequences
arg=sys.argv

#link the read to the barcode
barcodeReadFileName=arg[1]
seqOutFileName=arg[2]
seqIdFilterPct=float(arg[3])
#seqNamesFileName=arg[4]


#store seqNamesDict to allow less memory usage in following code
#seqNamesDict=dict()
#i=0
#with open(seqNamesFileName) as seqNamesFile:
#    for line in seqNamesFile:
#        line=line.split()
#        seqNamesDict[str(i)]=line[0]
#        i+=1


seqDict=dict()

totalNumBarcodes=0
totalNumDupBarcodes=0
totalNumBarcodesFilteredByPct=0

with open(barcodeReadFileName, "r") as barcodeReadFile:
    for line in barcodeReadFile:
        line=line.split()
        dupInfo=int(line[len(line)-1])
        totalNumBarcodes+=1
        
        if dupInfo>0:
            totalNumDupBarcodes+=1
            continue
        
        barcode=line[0]
        seqName=line[6]
        #seqName=seqNamesDict[seqId]
        
        barcodeArrStats=[ int(line[ind]) for ind in range(1, 6) ]
        
        totalReadLen=barcodeArrStats[1]
        numMismatches=barcodeArrStats[2]
        numGaps=barcodeArrStats[4]
        
        if float(numGaps+numMismatches)/totalReadLen > 1-seqIdFilterPct:
            totalNumBarcodesFilteredByPct+=1
            continue
        
        if seqName not in seqDict:
            seqDict[seqName]=[barcodeArrStats, [barcode]]

        else:
            seqDict[seqName][0]=[ seqDict[seqName][0][ind]+barcodeArrStats[ind] for ind in range(len(seqDict[seqName][0])) ]
            if barcode not in seqDict[seqName][1]:
                seqDict[seqName][1].append(barcode)
            
#del         
    
with open(seqOutFileName, "w") as seqOutFile:
    for seqName, info in seqDict.items():
        seqOutFile.write(seqName+"\t")

        #write out statistics
        for ind in range(len(info[0])):
            seqOutFile.write(str(info[0][ind])+"\t")

        #write out number of unique barcodes
        seqOutFile.write(str(len(info[1]))+"\t")

        #write out actual barcodes
        seqStr=""
        for ind in range(len(info[1])-1):
            seqStr+=str(info[1][ind])+":"
        seqStr+=str(info[1][len(info[1])-1])
        seqOutFile.write(seqStr+"\n")

with open(seqOutFileName+".stats", "w") as seqStatsOutFile:
    seqStatsOutFile.write(str(totalNumBarcodes)+"\t"+str(totalNumDupBarcodes)+"\t"+str(totalNumBarcodesFilteredByPct)+"\n")
