#!/bin/bash -v

#load modules samtool and minimap
module load samtools
module load minimap2/2.24

#minimap alignment to sam NP4-4
minimap2 -ax map-ont /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq.gz > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat.sam

#Converts back into bam, sorts and indexes
samtools view -bS /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat.sam > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat.bam

samtools sort -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat.bam

samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam

cd /share/lasallelab/Oran/nanopolish/f5c-v1.3

#f5c index with new fastq and more efficient caller
./f5c_x86_64_linux index -t 50 --iop 24 -d /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/fast5/NP4-4/ /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq

#cleans up any previous passed variants vcf to avoid downstream errors
cd /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/
rm /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/*.*
rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/tmp
rm -r /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/log


#activates clair3 env and loads samtools
conda activate clair3-1.0.4
module load samtools

#NP4-4 variant calling fastqconcats new bam from minimap using hg19 ref
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/hg19.fa --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3 --threads=20 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#filters for passed quality
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/NP4-4-PassedVariants.vcf

#indexes vcf file in prep for whatshap with bgzip then tabix
conda deactivate
conda activate /share/lasallelab/Oran/dovetail/luhmes/merged/oj
bgzip -i /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/NP4-4-PassedVariants.vcf
tabix -p vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/NP4-4-PassedVariants.vcf.gz

#activates whatshap env
conda deactivate
conda activate /share/lasallelab/Oran/miniconda3/whatshap-env

#phasing with whatshap (reference must be uncompressed and indexed)
whatshap phase --ignore-read-groups --reference /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/NP4-4-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/clair3/NP4-4-PassedVariants.vcf.gz /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam

#cleans up any previous methylation runs
rm /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-MethylationCall.tsv

#methylation calling with f5c (more efficient program than nanopolish *kept encountaring error so switched back to nanopolish)

#cd /share/lasallelab/Oran/nanopolish/f5c-v1.3

#./f5c_x86_64_linux call-methylation -t 60 --iop 30 -r /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq -b /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-MethylationCall.tsv

#Methylation calling using nanopolish
module load nanopolish
nanopolish call-methylation -t 30 -r /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-cat.fastq -b /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4x-cat_sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/hg19.fa > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-4-MethylationCall.tsv

# Prints this scary message after the ghost in the shell finishes running.  RIP: Zelda Rubinstein & Heather O'Rourke
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
