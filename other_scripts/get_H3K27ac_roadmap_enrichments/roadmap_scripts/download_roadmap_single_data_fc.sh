#!/bin/sh

roadmapID=$1
modType=$2
outDir=$3

wget "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/${roadmapID}-${modType}.fc.signal.bigwig" -O ${outDir}/${roadmapID}-${modType}.fc.signal.bigwig
