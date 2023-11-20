#!/bin/bash


#Nanoplot Nanopore QC plots


module load nanoplot/1.41.0
source activate nanoplot-1.41.0

#NP4-3
NanoPlot -t 10 --fastq /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.fastq.gz --tsv_stats --format pdf --legacy kde hex dot -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-summary-plots

#NP4-4
NanoPlot -t 10 --fastq /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq.gz --tsv_stats --format pdf --legacy kde hex dot -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat-summary-plots

#UDP4-2
NanoPlot -t 5 --fastq /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-2-cat.fastq.gz --tsv_stats --format pdf --legacy kde hex dot -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-2-cat-summary-plots

#UDP4-3
NanoPlot -t 5 --fastq /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-3-cat.fastq.gz --tsv_stats --format pdf --legacy kde hex dot -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-3-cat-summary-plots


#minimap alignment to sam

module load minimap2/2.24

#NP4-3
minimap2 -ax map-ont /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.sam

#NP4-4
minimap2 -ax map-ont /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.sam

#UDP4-2
minimap2 -ax map-ont /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-2-cat.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-2-cat.sam

#UDP4-3
minimap2 -ax map-ont /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-3-cat.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/UDP4-3-cat.sam


#sam to bam
module load samtools/1.15.1

samtools view -Sb  /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.sam > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.bam
