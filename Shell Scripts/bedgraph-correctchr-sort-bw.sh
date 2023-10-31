#!/bin/bash

#Changes directory to where binaries can be run
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed

#Fixes Chromosome names to 	UCSC *_CC*
./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/bedgraph/m_positive.bedgraph -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/bedgraph/NP4-3-m_positive_CC.bedgraph -a hg38.p13.chromAlias.txt

#removes 5th column from bedgraph file *_N5*
awk 'BEGIN { FS = OFS = "\t" } { $5 = ""; sub(/\t+$/, ""); print }' /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5.bedgraph

#liftover to hg19 *_hg19*
./liftOver /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5.bedgraph hg38ToHg19.over.chain.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph unMapped


#Sort and then merge overlapping regions(mean of score or max, sum, etc) *_SM*
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph.tmp

bedtools merge -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph.tmp -c 4 -d 0 -o max > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19_out.bedgraph

LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19_out.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19_SM.bedgraph


#changes to bigwig format *.bw*
./bedGraphToBigWig /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19_SM.bedgraph /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg19.chrom.sizes /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19_SM.bw


