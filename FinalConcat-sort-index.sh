#!/bin/bash

# Load samtools module
module load samtools

# Find BAM file with "concat" in the name in the current directory
bam_file=$(find . -name '*concat*.bam' -type f -print -quit)

# Check if a BAM file is found
if [ -z "$bam_file" ]; then
  echo "No BAM file with 'concat' in the name found."
  exit 1
fi

# Perform samtools sort on the BAM file
sorted_file="${bam_file%.bam}_sorted.bam"
samtools sort -o "$sorted_file" "$bam_file"

# Perform samtools index on the sorted file
samtools index "$sorted_file"


# Print the completed message after the code finishes running (It's huge so I can see the screen change from far away 8).  Bonus points if you get the reference!
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

TTTTTT H  H IIIII  SSS     H   H  OOO   U   U   SSS   EEEEE   IIIII  SSS      CCC  L    EEEE     A    N     N   !!
  TT   H  H   I   S    S   H   H O   O  U   U  S   S  E         I   S   S    C     L    E       A A   NN    N   !!
  TT   H  H   I    SS      H   H O   O  U   U   S     E         I    S       C     L    E      A   A  N N   N   !!
  TT   HHHH   I      S     HHHHH O   O  U   U    S    EEE       I      S     C     L    EEE   AAAAAAA N  N  N   !!
  TT   H  H   I       S    H   H O   O  U   U      S  E         I       S    C     L    E     A     A N   N N   !!
  TT   H  H   I   S   S    H   H O   O  U   U  S   S  E         I   S   S    C     L    E     A     A N    NN   
  TT   H  H IIIII  SSS     H   H  OOO    UUU    SSS   EEEEE   IIIII  SSS      CCC  LLLL EEEE  A     A N     N   !!
"