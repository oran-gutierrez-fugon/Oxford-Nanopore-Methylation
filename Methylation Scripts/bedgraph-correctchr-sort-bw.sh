#!/bin/bash

#may not be needed since have to sort after chromtoucsc with bedsort
bin$ ./sort-bed /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTEST.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsorted.bedgraph

bedToBigBed$ ./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsorted.bedgraph -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChr.bedgraph -a hg38.p13.chromAlias.txt

#Gave issues
./bedSort /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChr.bedgraph /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTbedsorted.bedgraph 


bedtools sort -faidx /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg38.p13.chrom.sizes -i m_posTESTsortedCorrectChr.bedgraph > m_posTESTsortedCorrectChrbedtoolssort.bedgraph


#I think was necessary but not 100% will try without previous bedtools sort and viceversa
LC_COLLATE=C sort -k1,1 -k2,2n /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssort.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATE.bedgraph



bedToBigBed$ ./bedGraphToBigWig /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATEno5thcol.bedgraph /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg38.p13.chrom.sizes /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATEno5thcol.bw


Expecting 4 words line 1 of /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATE.bedgraph got 5

#removes 5th row from bedgraph file
awk 'BEGIN { FS = OFS = "\t" } { $5 = ""; sub(/\t+$/, ""); print }' /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATE.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/m_posTESTsortedCorrectChrbedtoolssortCOLLATEno5thcol.bedgraph
