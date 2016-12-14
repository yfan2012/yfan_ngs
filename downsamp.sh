#!/bin/bash

#Downsample reads from a zipped fastq file

##List of files
data=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C*.fastq.gz

#Destination directory
destdir=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample50X/

#File name modifier for downsampled fastq
modif=down_50_

#number of reads you want to have per file at the end
num=3000000

for i in ${data}
do
    #Get rid of path stuff
    file=${i##*/}

    #Get rid of gz file extension
    name=${file%.*}

    #Downsample
    #-s is seed. Must be same for read pairs. I think it can be the same for all. 
    seqtk sample -s100 ${i} ${num} > ${destdir}/${modif}${name}
done





if [ 1 -eq 0 ]; then
    
r1776r1=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C1776_S1_L001_R1_001.fastq.gz
r1776r2=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C1776_S1_L001_R2_001.fastq.gz
r1939r1=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C1939_S2_L001_R1_001.fastq.gz
r1939r2=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C1939_S2_L001_R2_001.fastq.gz
r2898r1=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C2898_S3_L001_R1_001.fastq.gz
r2898r2=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C2898_S3_L001_R2_001.fastq.gz

seqtk sample -s100 ${r1776r1} 3000000 > ./down_50_1776r1.fastq &
seqtk sample -s100 ${r1776r2} 3000000 > ./down_50_1776r2.fastq &
seqtk sample -s105 ${r1939r1} 3000000 > ./down_50_1939r1.fastq &
seqtk sample -s105 ${r1939r2} 3000000 > ./down_50_1939r2.fastq &
seqtk sample -s110 ${r2898r1} 3000000 > ./down_50_2898r1.fastq &
seqtk sample -s110 ${r2898r2} 3000000 > ./down_50_2898r2.fastq

fi
