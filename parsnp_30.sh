#!/bin/bash

consensus=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/mpileup/fastas/
spades=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/spades/
ref=/mithril/Data/NGS/Reference/cneo/grubii/grubii.fasta
out=~/Dropbox/Lab/fungus_zhang/

~/software/parsnp/parsnp -r $ref -d $consensus -p 12 -o ${out}parsnp_30_cons
#~/software/parsnp/parsnp -r $ref -d $spades -p 12 -o ${out}parsnp_30_spades

mv ${out}parsnp_30_cons/parsnp.ggr ${out}parsnp_30_cons/parsnp_cons_30.ggr
#mv ${out}parsnp_30_spades/parsnp.ggr ${out}parsnp_30_spades/parsnp_spades.ggr 
