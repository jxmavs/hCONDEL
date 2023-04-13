#!/bin/sh

wiggleFileName=$1
outFile=$2

bigWigToBedGraph ${wiggleFileName} ${outFile}
#wget "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/${roadmapID}-${modType}.pval.signal.bigwig" -O ${outDir}/${roadmapID}-${modType}.pval.signal.bigwig
