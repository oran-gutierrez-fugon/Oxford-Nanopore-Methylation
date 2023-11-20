#!/bin/bash -v

#Author: Oran Gutierrez Fugon MD PhD, LaSalle Lab, Segal Lab, Integrative Genetics and Genomics graduate group UC Davis

#Although this is structured as a shell script I would recommend each section to be run individually to deal with errors as they arise. Also ran into issues switching between conda env while using screen and some steps are only included to clean up previous attempts and may not be necessary.

#Generally have not seen any part of this pipeline taking up more than 12 GB of memory with 60 cores going at a time but to be safe and respect epigenerate use the this command to limit ram usuage before killing the job at 75GB. Core usage options will vary with resorces available on epigenerate at the time of running.
ulimit -v 75000000

#load modules samtool and minimap
module load samtools
module load minimap2/2.24

#concatenates all basecalled fastqs from all flushes to catcat folder (do not use methylation fastqs) see nanomethphase github
#May want to remove individual flush concat fastq at the end but recommend leave them until you finish the pipeline in case you need to process flushes individualy and then merge the methylation call tsv
#*****WARNING***** input fastq preflush,PF,PF2 folders must be verified with cd one by one not just find/replace replicate name
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05182023/20230518_1823_3G_PAO36704_d713bcfc/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-.fastq.gz 
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05192023_PF/20230519_1725_3G_PAO36704_22df49af/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF.fastq.gz 
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05202023_PF2/20230520_1804_3G_PAO36704_e7bf9f77/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF2.fastq.gz

cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-.fastq.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF.fastq.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF2.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz

#converts fast5 raw files to more efficient blow5 and puts all flushes in same directory. slow5tools was installed in base conda env
#*****WARNING***** input fast5 preflush,PF,PF2 folders must be verified with cd one by one not just find/replace replicate name
conda activate
slow5tools fast5toslow5 -p 60 --to blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05182023/20230518_1823_3G_PAO36704_d713bcfc/fast5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05192023_PF/20230519_1725_3G_PAO36704_22df49af/fast5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05202023_PF2/20230520_1804_3G_PAO36704_e7bf9f77/fast5 -d /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3; echo "It is done" | mail -s "UDP4-3 fast5toslow5" ojg333@gmail.com

#merges slow5 dir to 1 single blow5 file
mkdir /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/
slow5tools merge -t 60 --to blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3  -o /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5

#Cleans up single blow5 files since they are no longer needed
rm /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/*.slow5 

# activates f5c env and indexes blow5 and fastq
conda activate f5c
f5c index -t 60 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz

#minimap2 aligment + samtools sam to bam + indexes bam
minimap2 -a -x map-ont -t 60 -2 /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz | samtools sort -T tmp -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam
samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam; echo "It is done.." | mail -s "UDP4-3 samtools index done" ojg333@gmail.com

#cleans up any previous methylation runs if needed
#rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv

#methylation calling with f5c (more efficient program than nanopolish) must include option --pore r10 for r10 chemistry (thanks Logan!)
f5c call-methylation --pore r10 -x hpc-high --meth-out-version 2 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5 -r /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz  -b /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv; echo "It is done.." | mail -s "UDP4-3 f5c call-methylation done" ojg333@gmail.com

#cleans up any previous passed variants vcf to avoid mixed data errors
#cd /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/
#rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/*.*
#rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/tmp
#rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/log

#activates clair3 env and loads samtools
conda activate clair3-1.0.4
module load samtools

#UDP4-3 variant calling fastqconcats new bam from minimap using hg19 ref. Model is dependent on chemistry and basecaller used (see clair3 and reiro githubs)
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/hg19.fa --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3 --threads=60 --platform=ont --model_path=/share/lasallelab/Oran/dovetail/luhmes/methylation/clair3model/rerio/clair3_models/r1041_e82_400bps_sup_g615

#filters for passed quality leave unzipped for now since will need to be in bgzip compression format for downstream steps
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf

#indexes vcf file in prep for whatshap with bgzip then indexes with tabix. Using HiChIP oj conda env with bgzip and tabix already installed 
conda deactivate
conda activate /share/lasallelab/Oran/dovetail/luhmes/merged/oj
bgzip -i /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf
tabix -p vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf.gz

#switches to whatshap conda env
conda deactivate
conda activate /share/lasallelab/Oran/miniconda3/whatshap-env

#phasing with whatshap (reference must be uncompressed and indexed).  Since whatshap is not optimized for multicore this is a good step to do in parallel with other samples on the cluster after running clair3 separately or together with clair3 running for the next sample. May also substitute illumina variant vcf file if available.
whatshap phase --ignore-read-groups --reference /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam

#activates nanomethphase conda env and python step 1 methyl call processor + index. Must be in Nanomethphase folder (Thanks Osman!)
conda deactivate
conda activate /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase/nanometh-environment
cd /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase
python nanomethphase.py methyl_call_processor -mc /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv -t 60 | sort -k1,1 -k2,2n -k3,3n | bgzip > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.bed.gz && tabix -p bed /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.bed.gz

#phases the methylome (the moment you've all been waiting for)
python nanomethphase.py phase --include_indels -b /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam -v /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-whatshap_phased.vcf -mc /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.bed.gz -r /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome -of bam,methylcall,bam2bis -t 50

#aggregates data from both strands (requieres datamash installation)
#use correct file names from previous step, differential methylation in next step does this automatically so can skip
conda activate datamash
sed '1d' /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP1_MethylFrequency.tsv | awk -F'\t' '{if ($4=="-") {$2=$2-1;$3=$3-1}; print $1,$2,$3,$5,$6}' OFS='\t' | sort -k1,1 -k2,2n | datamash -g1,2,3 sum 4,5 | awk -F'\t' '{print $0,$5/$4}' OFS='\t' | sed '1i chromosome\tstart\tend\tNumOfAllCalls\tNumOfModCalls\tMethylFreq' > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/UDP4-3_aggregated_HP1_MethylFrequency.tsv

sed '1d' /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP2_MethylFrequency.tsv | awk -F'\t' '{if ($4=="-") {$2=$2-1;$3=$3-1}; print $1,$2,$3,$5,$6}' OFS='\t' | sort -k1,1 -k2,2n | datamash -g1,2,3 sum 4,5 | awk -F'\t' '{print $0,$5/$4}' OFS='\t' | sed '1i chromosome\tstart\tend\tNumOfAllCalls\tNumOfModCalls\tMethylFreq' > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/UDP4-3_aggregated_HP2_MethylFrequency.tsv

#creates directory and cleans up previous dma tries
mkdir /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA
rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/*.*

#If running script linearly will need to go back to nanomethphase conda env
conda deactivate
conda activate /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase/nanometh-environment
cd /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase

#Differential methylation analysis (if you've made it this far, lets go a little farther)
#Check folders and file names match with previous steps but not datamash output since this will aggregate automatically
#see DSS ddocumentation for all options and output file format
#Had to install sys for commandline R in nanomethphase env using R then install.packages("sys")
python nanomethphase.py dma -c 1,2,4,5,7 -ca /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP1_MethylFrequency.tsv -co /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP2_MethylFrequency.tsv -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/ -op DMA

#Visualization for viewing in UCSC genome browser:

#Converts output tsv files to bedgraph 4 column format (can take the read count column instead of methylation if want to make a coverage plot track)
awk 'BEGIN {FS="\t"; OFS="\t"}
NR > 1 {print $1, $2, $3, $7*100}' /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP1_MethylFrequency.tsv > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP1_MethylFrequency.bedGraph

awk 'BEGIN {FS="\t"; OFS="\t"}
NR > 1 {print $1, $2, $3, $7*100}' /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_methylome_NanoMethPhase_HP2_MethylFrequency.tsv > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP2_MethylFrequency.bedGraph

#For combining replicates from begraphs (Not recommended since can show both replicates overlayed using trackhubs)
#cat replicate1.bedGraph replicate2.bedGraph | sort -k1,1 -k2,2n > combined_sorted.bedGraph

#Changes directory to where bedGraphToBigWig binaries can be run
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed

#sorts bedgraph with bedSort
./bedSort /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP1_MethylFrequency.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP1_MethylFrequency_sorted.bedGraph

./bedSort /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP2_MethylFrequency.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP2_MethylFrequency_sorted.bedGraph

#Changes to bigwig format (.bw) for fast viewing in UCSC genome browser
./bedGraphToBigWig /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP1_MethylFrequency_sorted.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg19.chrom.sizes /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP1_MethylFrequency.bw

./bedGraphToBigWig /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP2_MethylFrequency_sorted.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg19.chrom.sizes /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/bedgraphs/UDP4-3_HP2_MethylFrequency.bw


#Visualization of Differential Methylation Analysis (DMA)
#Convert space delimited txt callDMR file to tab delimited bedGraph with column 8 converted to percent by multiplying times 100

rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR.bedGraph
rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR_Sorted.bedGraph

#Converts output DMA txt files to bedgraph 4 column format (can take the read count column instead of methylation if want to make a coverage plot track)
awk 'BEGIN {FS="\t"; OFS="\t"}
NR > 1 {print $1, $2, $3, $8*100}' /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/DMA_callDMR.txt > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR.bedGraph


cd /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed

#sorts bedgraph with bedSort
./bedSort /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR_Sorted.bedGraph


#Changes to bigwig format (.bw) for fast viewing in UCSC genome browser
./bedGraphToBigWig /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR_Sorted.bedGraph /share/lasallelab/Oran/dovetail/luhmes/methylation/bedToBigBed/hg19.chrom.sizes /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/DMA/UDP4-3_callDMR_Sorted.bw


#To view bw files just upload to bioshare, copy link, and create track hub on UCSC genome browser 

# Prints this scary message after the ghost in the shell finishes running so my lazy bones can see it finish from far away.  Bonus points if you can figure out the reference, RIP: Zelda Rubinstein & Heather O'Rourke
echo "
                                                     @@@@@@@@
                                             %&&&&@&          @@@@@(@&@&&#                                                             
                                         /@&%                            .&&@.                                                       
                                      *&@                                      #&&                                                    
                                   .&&                                            /@%                                                 
                                 %&                                                 .&@                                               
                               #&                                                    ..&&                                             
                              &.                                                       ,/&                                            
                            #&                                                          ,.&/                                          
                           @#                                                            ..&#                                         
                          @/                                                              ,.&*                                        
                         &(                                                               .,.@                                        
                        %&                                                                 ,.&&                                       
                        &                                                                  ...&                                       
                       @           @(&&&&&#                @&&&&&@*                        .,.&&                                      
                      (&         ,%&&&&&&&&&@           (&&&&&&&&&&,,                      .,,(&                                      
                      &         .*&&@#&&&&&&&@         #&&&&@&&&&&&&,,                     .,,,&                                      
   .&&&@@@&&@#@&@/   #%         .*&&&&&&&&&&&&         &&&&&&&&&&&&&%,                     .,,.&    #&&@#&&&@@@#&&%.                   
   &.             &#&&           .&&&&&&&&&&&          *&&&&&&&&&&&&,,                     ,,,.&&###         .,,. %@                  
    @% ,            ,&            *.#&&&&&&              &&&&&&&&@&,.                      ,&,             .,,,. %&                   
      &@ ,.        .,                                       &&&                                          ,,,,  @&                     
        &% .,      ,.                                                                                  ,,,,  @@                       
          @* ... .,,,             ,###,       /&&&&&&\      ,###,                                   .,,,,  @&                         
           #& ,,,,,,.          @&&&&&&&&&&&&&&&&&&&&&&&&&&@&&&&&&&&&&&&                           ..,,,. @&                           
             &  .,,.          &%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&,%&                        ,,,,, *&                             
              &# ,,(          @#&&&&&&&&&&&&@&&&&&**&&&@&&&&&&&&&&&&&&&.&                      ,,,,,. @%                             
               && ,&            &&&,   '&&&&'          '&&&&'     &&&&#                       ,,.,...&                             
                %@ &                                                                        .*,,,,.,@                     
                 &%&                                                                        ,,,,,.(@                              
                  @&                                                                        ,,..,*&                       
                  ,&                                                                     ..,.,., &&                 
                  /@                                                              ,        ,,,,,@ /@                                  
                  %&                          ,                                   .        ,,,,.,, &(                                 
                  @(                                                              #.       ,,,,,,,, @.                                
                  &.                         .                                    @,.      ,.   .,,. &                                
                  &                          &.                                   &,,            .,,. &                               
                 &@                          #.           .                       *%.,            ,,,  &                              
                 @/                          *,           ,                        &.,.            .,.  &                             
                .&            &             /.,           ,.                       @/,.             ,,,  &                            
               &@           .&             &..           .@                        &,,,            .,,,  @.                          
               .@           ,@              &..          .,&                        @*,,,            ,,,,. &*                         
               &&          ..@              &,.          ,.&                         &,,,,            ,,,,. @.                        
              /&          .,&              .@.,          ,,#.                        &#.,,,            ,,,,, &                        
              @,         .,,@              /@.,          ,,*#                         &.,,,,            ,,,.&%                        
             @@         .,,&               %%,,,         ,,.@                         ,@,,,,..           ,#&                          
            (&         .,,%&               &(,,,        .,,,&                          &%,,,,,,,       @#,                            
            &,        ,,..&                @*,,,        .,,,@                           &/.,.,,.,.%@@.                                
            &&      ,,,,,@*                &,,,.        ,,,.%(                           &,,,,&@/                                     
              #%&#.,,,,,*&                 &,,,,,       ,,,.,&                            &@&                                         
                    &&%,&,                 &*,.,.       ,,,,,&                            .&                                          
                       %@                  &(,,,...   ..,,,,,%&                            *&                                         
                       &%                  %@,    '''''      .&                             &&                                        
                       &&                   &                 &&                            .@                                        
                         &@                @@                  @*                         /&*                                         
                           %%%%&&&&&&&&@&@*                      &@                     &&%                                            
                                                                   (@&&(&&&&&&&&&&&&(%@&&*  

TTTTTT H  H IIIII  SSS     H   H  OOO   U   U   SSS   EEEEE   IIIII  SSS      CCC  L    EEEE     A    N     N
  TT   H  H   I   S    S   H   H O   O  U   U  S   S  E         I   S   S    C     L    E       A A   NN    N
  TT   H  H   I    SS      H   H O   O  U   U   S     E         I    S       C     L    E      A   A  N N   N
  TT   HHHH   I      S     HHHHH O   O  U   U    S    EEE       I      S     C     L    EEE   AAAAAAA N  N  N
  TT   H  H   I       S    H   H O   O  U   U      S  E         I       S    C     L    E     A     A N   N N
  TT   H  H   I   S   S    H   H O   O  U   U  S   S  E         I   S   S    C     L    E     A     A N    NN
  TT   H  H IIIII  SSS     H   H  OOO    UUU    SSS   EEEEE   IIIII  SSS      CCC  LLLL EEEE  A     A N     N
"
