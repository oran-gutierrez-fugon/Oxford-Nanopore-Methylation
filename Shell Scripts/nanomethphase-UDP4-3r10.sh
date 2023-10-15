#!/bin/bash -v

#Author: Oran Gutierrez Fugon MD PhD, LaSalle Lab, Segal Lab, Integrative Genetics and Genomics graduate group UC Davis

#Although this is structured as a shell script I would recommend each section to be run individually to deal with errors as they arise

#load modules samtool and minimap
module load samtools
module load minimap2/2.24

#concatenates all basecalled fastqs from all flushes to catcat folder (do not use methylation fastqs) see nanomethphase github
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05182023/20230518_1823_3G_PAO36704_d713bcfc/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-.fastq.gz 
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05192023_PF/20230519_1725_3G_PAO36704_22df49af/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF.fastq.gz 
cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05202023_PF2/20230520_1804_3G_PAO36704_e7bf9f77/fast5/basecalling/pass/*.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF2.fastq.gz

cat /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-.fastq.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF.fastq.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3-PF2.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz

#May want to remove individual flush concat fastq but recommend leave them until you finish the pipeline in case you need to process flushes individualy and then merge the tsv

#converts fast5 raw files to more efficient blow5 and puts all flushes in same directory. slow5tools was installed in base conda env
#*****WARNING***** input fast5 preflush,PF,PF2 folders must be verified with cd one by one not just find/replace replicate name
conda activate
slow5tools fast5toslow5 -p 60 --to blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05182023/20230518_1823_3G_PAO36704_d713bcfc/fast5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05192023_PF/20230519_1725_3G_PAO36704_22df49af/fast5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/PROM0151_LaSalle_NP4_4_05202023_PF2/20230520_1804_3G_PAO36704_e7bf9f77/fast5 -d /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3; echo "It is done" | mail -s "UDP4-3 fast5toslow5" ojg333@gmail.com

#merges slow5 dir to 1 single blow5 file
mkdir /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/
slow5tools merge -t 60 --to blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3  -o /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5

#Cleans up single blow5 files since they are no longer needed
#rm /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/*.slow5 

# activates f5c env and indexes blow5 and fastq
conda activate f5c
f5c index -t 60 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz

#minimap2 aligment + samtools sam to bam + indexes bam
minimap2 -a -x map-ont -t 60 -2 /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa.gz /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz | samtools sort -T tmp -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam
samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam; echo "UDP4-3 samtools index done title" | mail -s "UDP4-3 samtools index done content" ojg333@gmail.com

#cleans up any previous methylation runs if needed
#rm /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv

#methylation calling with f5c (more efficient program than nanopolish) must include option --pore r10 for r10 chemistry (thanks Logan!)
f5c call-methylation --pore r10 -x hpc-high --meth-out-version 2 --slow5 /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/slow5/UDP4-3/blow/UDP4-3cat.blow5 -r /share/lasallelab/Oran/dovetail/luhmes/nanoRAW/catcat/UDP4-3.fastq.gz  -b /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv; echo "UDP4-3 f5c call-methylation done" | mail -s "UDP4-3 f5c call-methylation done" ojg333@gmail.com

#cleans up any previous passed variants vcf to avoid downstream errors
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-3/clair3/
rm /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-3/clair3/*.*
rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-3/clair3/tmp
rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/UDP4-3/clair3/log

#activates clair3 env and loads samtools
conda activate clair3-1.0.4
module load samtools

#UDP4-3 variant calling fastqconcats new bam from minimap using hg19 ref
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/hg19.fa --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3 --threads=60 --platform=ont --model_path=/share/lasallelab/Oran/dovetail/luhmes/methylation/clair3model/rerio/clair3_models/r1041_e82_400bps_sup_g615

#filters for passed quality leave unzipped for now since will need to be in bgzip compression format for downstream steps
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf

#indexes vcf file in prep for whatshap with bgzip then indexes with tabix
conda deactivate
conda activate /share/lasallelab/Oran/dovetail/luhmes/merged/oj
bgzip -i /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf
tabix -p vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf.gz

#activates whatshap env
conda deactivate
conda activate /share/lasallelab/Oran/miniconda3/whatshap-env

#phasing with whatshap (reference must be uncompressed and indexed)
whatshap phase --ignore-read-groups --reference /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa -o /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3/clair3/UDP4-3-PassedVariants.vcf.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3_sorted.bam

#activates nanomethphase conda env and python step 1 methyl call processor (Thanks Osman!)
conda deactivate
conda activate /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase/nanometh-environment
cd /share/lasallelab/Oran/test_nanomethphase/NanoMethPhase
python nanomethphase.py methyl_call_processor -mc /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.bed.gz && tabix -p bed /share/lasallelab/Oran/dovetail/luhmes/methylation/phasing/UDP4-3-MethylationCall.bed.gz


# Prints this scary message after the ghost in the shell finishes running.  Bonus points if you get the reference, RIP: Zelda Rubinstein & Heather O'Rourke
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
