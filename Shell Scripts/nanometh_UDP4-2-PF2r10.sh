#!/bin/bash -v

# Test UDP4-2-PF2 (non multirun cat)
module load samtools
module load minimap2/2.24

# changes folder & cat fastq (specific folders vary)
cd /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_UDP4-2_05042023_PF/20230504_1405_2G_PAM05619_b1291f1b/fast5/basecalling/pass/
cat *.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-2-PF.fastq.gz

# create slow5 (or blow5) from single read fast5 (slow5tools was installed in my base miniconda env)
conda activate
slow5tools fast5toslow5 -p 60 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_UDP4-2_05052023_PF2/20230505_1817_2G_PAM05619_b73e9ab9/fast5 -d /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/PF2

# must merge multiple slow5 (or binary blow5)
slow5tools merge -t 60 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/PF2 -o /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/PF2/UDP4-PF2.blow5

# remove multiple slow5 (or blow5)
cd /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/PF2/
rm *.slow5

# index using optimized fork of nanopolish f5c that supports R10 chemistry
conda activate f5c
f5c index -t 60 --pore r10 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/blow/UDP4-2cat.blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-2.fastq.gz 

# maps to ref creates sorts and indexes bam
minimap2 -a -x map-ont /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-2.fastq.gz | samtools sort -T tmp -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2_sorted.bam
samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2_sorted.bam

#f5c more efficient caller, use slow5 directories and concatenated fastq used in minimap step
f5c call-methylation -t 50 --pore r10 -x hpc-high --meth-out-version 2 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-2/PF2/UDP4-PF2.blow5 -r /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-2-PF2.fastq.gz  -b /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2-PF2_sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2-PF2-MethylationCall.tsv

# nanomethphase conda env and python step 1 methyl call processor
conda deactivate
conda activate /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase/nanometh-environment
cd /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase
python nanomethphase.py methyl_call_processor -mc /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2-PF2-MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2-PF2-MethylationCall.bed.gz && tabix -p bed /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-2-PF2-MethylationCall.bed.gz

