# Oxford Nanopore phased methylation in LUHMES neurons

Allele-specific (phased) CpG methylation from Oxford Nanopore long-read whole-genome sequencing of LUHMES cells, comparing undifferentiated precursors (UDP4) to day-7 neurons (NP4). The goal was to find differentially methylated regions at the 15q11-q13 locus that set the developmental boundary controlling *UBE3A-ATS*, resolved by parental allele.

Part of the methods behind:

> Gutierrez Fugón OJ, Sharifi O, Heath NC, et al. *Integration of CTCF Loops, Methylome, and Transcriptome in Differentiating LUHMES as a Model for Imprinting Dynamics of the 15q11-q13 Locus in Human Neurons.* Hum Mol Genet. 2024. https://doi.org/10.1093/hmg/ddae111

## Genome browser sessions

Interactive UCSC Genome Browser sessions with the CTCF loops, Oxford Nanopore phased methylation, and RNA-seq tracks for the 15q11-q13 region in LUHMES:

- **AS imprinted region (Figure 7A)**: CTCF loops, Oxford nanopore phased methylation and RNAseq in LUHMES cell line for the 15q11-q13 region. View of AS imprinted region as shown in Figure 7A of Gutierrez Fugon et al, Human Molecular Genetics, 2024. https://genome.ucsc.edu/s/ojg333/AS_Locus_Fig7A_Fugon_2024
- **MAGEL2 loop anchor, zoomed (Figure 9A)**: CTCF loops, Oxford nanopore phased methylation and RNAseq in LUHMES cell line for the 15q11-q13 region. Zoomed view of the MAGEL2 neuron specific loop anchor as shown in Figure 9A of Gutierrez Fugon et al, Human Molecular Genetics, 2024. https://genome.ucsc.edu/s/ojg333/MAGEL2_Fig9A%20_Fugon_2024

## Main result

A paternally hypomethylated DMR near the *SNRPN* anchor appears only in neurons, while paternally hypermethylated DMRs near *PWAR1* appear only in undifferentiated cells. These allele-specific patterns line up with the neuron-specific CTCF loops and mark where *UBE3A-ATS* transcription is allowed to extend.

## Pipeline

Sample labels: **NP4** = neurons, **UDP4** = undifferentiated. The phasing scripts live in `Shell Scripts/`, with earlier iterations kept in `Shell Scripts/OLD/`.

**1. Merge per-flush BAMs** (`Concatenate all Flushes.sh`, `ConcatenateGuppyBams-Pass.sh`, `FinalConcat-sort-index.sh`)
Each Nanopore flush produces its own pass BAMs. They are merged header-aware so the methylation tags survive, then `samtools sort` and `index`.

**2. Align and call variants** (`Shell Scripts/OLD/minimapFASTQtoBAM*.sh`, `Shell Scripts/OLD/clair3*.sh`)
Reads are aligned with minimap2 and variants called with Clair3 to provide the heterozygous sites used for phasing.

**3. Phase methylation by allele** (`Shell Scripts/nanomethphase-*.sh`, `Shell Scripts/OLD/methylationphasing-*.sh`)
nanomethphase (with whatshap and f5c) assigns reads and their CpG calls to haplotype 1 or 2, giving parent-of-origin methylation for each sample (NP4-3, NP4-4, UDP4-2, UDP4-3).

**4. Pile up methylation** (`modkit-pileup.sh`, `Shell Scripts/modkit-pileupBEDGRAPH.sh`)
`modkit pileup` produces per-CpG methylation bedGraphs.

**5. Fix chromosome names and lift over** (`ChrCorrector.sh`, `correctchromnames.r`, `hg38ToHg19.over.chain.gz`, `Shell Scripts/bedgraph-correctchr-sort-bw.sh`)
A common quiet failure is the bedGraph and reference disagreeing on chromosome names, which silently drops rows downstream. `ChrCorrector.sh` reads a name-mapping table into an associative array and rewrites the first column. The liftover chain and `bedgraph-correctchr-sort-bw.sh` then move tracks to hg19, sort, and convert to bigWig for the browser.

**6. DMR follow-up** (`Shell Scripts/GOterms*.r`)
GO-term enrichment on the differentially methylated regions, split by cell state and ontology (biological process for neurons and undifferentiated, cellular component, Reactome).

## Phased track hubs

`phased_trackhubs/` holds ready-to-load UCSC track hubs of the phased methylation, separated by haplotype: `Neurons-H1`, `Neurons-H2`, `Undifferentiated-H1`, `Undifferentiated-H2`, plus combined `DMA-Neurons` and `DMA_Undif` differential-methylation hubs. Each carries its own `hub.txt`, `genomes.txt`, and hg19 `trackDb.txt`.

## Data upload

`Shell Scripts/SRAupload/*split.sh` split the per-sample reads for submission to SRA.

## Notes

Scripts were written for a SLURM cluster, so `module load` lines and absolute paths need editing for another setup. Inputs are methylation-tagged BAMs from the Nanopore basecaller.
