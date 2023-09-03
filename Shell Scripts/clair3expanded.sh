#!/bin/bash

#Original from github pattern
./run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3 --threads=20 --platform=ont --model_path=/share/lasallelab/Oran/dovetail/luhmes/methylation/clair3model/r941_prom_sup_g5014

#better tuned and working but error 
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3 --threads=50 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#Demo working correctly
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/demo/HG003_chr20_demo.bam --ref_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/demo/GRCh38_no_alt_chr20.fa --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/demo/output --threads=50 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#Tried with NP4-3 bam and the demo ref put contig name chr1, chr2 etc error due to bam file names
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/demo --threads=50 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#to correct chromosome names, did not fix error, probably due to bam format
samtools view -h /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.bam | ./chromToUcsc -a hg38.p13.chromAlias.txt | samtools -bS > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat_cc.bam

#previously used chromToUcsc version
bedToBigBed$ ./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/bedgraph/m_positive.bedgraph -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/bedgraph/NP4-3-m_positive_CC.bedgraph -a hg38.p13.chromAlias.txt

#to convert to sam first WORKED
samtools view -h /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.bam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.sam

#chromToUcsc for sam WORKED
bedToBigBed$ ./chromToUcsc -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat.sam -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam -a hg38.p13.chromAlias.txt

#liftover to hg19 (failed:  Expecting integer field 2 line 1 of /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam, got VN:1.6)
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed
./liftOver /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam hg38ToHg19.over.chain.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC_hg19.sam unMapped
#did not do
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/correctedhg19/NP4-3-m_positive_CC_N5_hg19.bedgraph.tmp

#Turns back into bam
samtools view -bS /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.sam > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCC.bam

#sort chromosome corrected CC bam then index
samtools sort -o NP4-3concatCCsorted.bam NP4-3concatCC.bam
samtools index NP4-3concatCCsorted.bam

#concatenates sequencing summary file to enable and speed up indexing
cat /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05192023_PF/20230519_1725_3E_PAO28520_cbd61198/fast5/basecalling/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05202023_PF2/20230520_1804_3E_PAO28520_6c925a53/fast5/basecalling/sequencing_summary.txt > /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/seqsummarycat/sequencing_summary.txt

#index fastq files with fast5 for just preflush NP4-3 (redownload using old nanopolish) WORKED!
nanopolish index -d /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/ -s /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/pass/concat/NP4-3.fastq

#index fastq files with fast5 for concat NP4-3
nanopolish index -d /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/fast5/NP4-3/ -s /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/seqsummarycat/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/NP4-3-cat.fastq.gz





#methylation calling (failed since no index db file)
nanopolish call-methylation -t 40 -q cpg -r /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/pass/concat/NP4-3.fastq -b /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-MethylationCall.tsv

#NP4-3 chromosome corrected CC bam attempt to skip to step 2 with guppy generated bam after chromosome correction using another ref genome (now working but try with correct ref .p13 as well)
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect --threads=30 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#filters for passed Quality
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > PassedVariants.vcf

#phasing
whatshap phase --ignore-read-groups --indels --reference /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/NP4-3-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/chromcorrect/PassedVariants.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam






#NP4-3 chromosome corrected CC bam attempt to skip to step 2 with guppy generated bam after chromosome correction using another ref genome (now working but try with correct ref .p13 as well)
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concatCCsorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/refchromcorrect --threads=30 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014



#identical
diff -sq /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/pass/fastq_runid_a427361483a5cf4567d4e4a308f56521c8d0df8c_0_0.fastq.gz /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/pass/fastq_runid_a427361483a5cf4567d4e4a308f56521c8d0df8c_0_0.fastq.gz

#dif
diff -sq /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/redownload/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/5mc-5hmc/pass/fastq_runid_a427361483a5cf4567d4e4a308f56521c8d0df8c_0_0.fastq.gz 




@HD     VN:1.6  SO:coordinate


