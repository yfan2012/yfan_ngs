#!/bin/bash

#Count the number of reads in a zipped fastq file. 
fastqs=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C*

for i in $fastqs
do
    echo $i
    linecount=$(gunzip -c $i | wc -l)
    ##echo ${linecount}
    expr ${linecount} / 4
done
