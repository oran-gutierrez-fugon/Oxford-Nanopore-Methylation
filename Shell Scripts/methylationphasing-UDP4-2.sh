#!/bin/bash -v

conda activate clair3-1.0.4
module load samtools
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed

#Convert to concat bam to sam first
samtools view -h /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concat.bam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concat.sam

#chromToUcsc for sam corrects names from hg38.p13
./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concat.sam -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCC.sam -a hg38.p13.chromAlias.txt

#Turns back into bam
samtools view -bS /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCC.sam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCC.bam

#sort chromosome corrected CC bam then index
samtools sort -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCCsorted.bam /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCC.bam
samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCCsorted.bam

#UDP4-2 chromosome corrected CC bam variant calling using GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCCsorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/clair3 --threads=30 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#filters for passed Quality
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/clair3/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/clair3/PassedVariants.vcf

#change and activate whatshap environment
conda deactivate
conda activate /share/lasallelab/Oran/miniconda3/whatshap-env
module load samtools

#phasing
whatshap phase --ignore-read-groups --reference /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/clair3/UDP4-2-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/clair3/PassedVariants.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-2/UDP4-2concatCCsorted.bam