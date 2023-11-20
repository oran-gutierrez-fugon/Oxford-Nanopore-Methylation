#!/bin/bash

./run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/NP4-3concat_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_genomic.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3 --threads=20 --platform=ont --model_path=/share/lasallelab/Oran/dovetail/luhmes/methylation/clair3model/r941_prom_sup_g5014