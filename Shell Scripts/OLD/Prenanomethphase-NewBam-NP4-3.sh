#!/bin/bash -v

#concatenates sequencing summary file to enable and speed up indexing
cat /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05182023/20230518_1823_3E_PAO28520_86921a60/fast5/basecalling/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05192023_PF/20230519_1725_3E_PAO28520_cbd61198/fast5/basecalling/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/PROM0151_LaSalle_NP4_3_05202023_PF2/20230520_1804_3E_PAO28520_6c925a53/fast5/basecalling/sequencing_summary.txt > /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/seqsummarycat/sequencing_summary.txt

#activates nanopolish env and loads samtools
conda activate nanopolish
module load samtools

#index fastq files with fast5 for concat NP4-3
nanopolish index -d /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/fast5/NP4-3/ -s /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/seqsummarycat/sequencing_summary.txt /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/NP4-3-cat.fastq.gz

#sorts and indexes minimap created bam file from fastqconcats directory (already done for NP4-3)
samtools sort -o /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-sorted.bam /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat.bam
samtools index /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-sorted.bam

#methylation calling (failed since no index db file should work now that index may be working for concat, but should wait for new install)
nanopolish call-methylation -t 40 -q cpg -r /share/lasallelab/Oran/dovetail/luhmes/nanoporeRAW/NP4-3-concatpass/NP4-3-cat.fastq.gz -b /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-sorted.bam -g /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-MethylationCall.tsv

#activates clair3 env and loads samtools
conda deactivate
conda activate clair3-1.0.4
module load samtools

#NP4-3 variant calling fastqconcats bam from minimap using GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
run_clair3.sh --bam_fn=/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-sorted.bam --ref_fn=/share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --output=/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3 --threads=30 --platform=ont --model_path=/software/anaconda3/23.1.0/lssc0-linux/envs/clair3-1.0.4/bin/models/r941_prom_sup_g5014

#filters for passed quality
gunzip -c /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/merge_output.vcf.gz | awk '$1 ~ /^#/ || $7=="PASS"' > /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/PassedVariants.vcf

#activates whatshap env and loads samtools
conda deactivate
conda activate /share/lasallelab/Oran/miniconda3/whatshap-env
module load samtools

#phasing with whatshap
whatshap phase --ignore-read-groups --reference /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/NP4-3-whatshap_phased.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-3/clair3/PassedVariants.vcf /share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/NP4-3-cat-sorted.bam


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
