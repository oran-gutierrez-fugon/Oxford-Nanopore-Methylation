#!/bin/bash

#Convert to concat bam to sam first
samtools view -h /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.bam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.sam

#chromToUcsc for sam corrects names from hg38.p13
bedToBigBed$ ./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.sam -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam -a hg38.p13.chromAlias.txt

#Turns back into bam
samtools view -bS /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.bam

#sort chromosome corrected CC bam then index
samtools sort -o NP4-3concatCCsorted.bam NP4-3concatCC.bam
samtools index NP4-3concatCCsorted.bam

#NP4-3 chromosome corrected CC bam variant calling using GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect --threads=30 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#filters for passed Quality
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > PassedVariants.vcf

#phasing
whatshap phase --ignore-read-groups --indels --reference /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/NP4-3-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/PassedVariants.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam