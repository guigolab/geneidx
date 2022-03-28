#!/bin/bash

##########################
# Building coding statistic (Markov model order 5)
#Generating the initial and transition matrices
# without stop codon
#################
#markov5
##################

gawk '{print \$1,substr(\$2,1,length(\$2)-3)}' all.cds_filter1.tbl | gawk -f MarkovMatrices.awk 5 set1.cds

sort +1 -2  -o set1.cds.5.initial set1.cds.5.initial
sort +1 -2  -o set1.cds.5.transition set1.cds.5.transition

gawk -f MarkovMatrices-noframe.awk 5  set1.intron all.intron_filter1.tbl

sort -o set1.intron.5.initial set1.intron.5.initial
sort -o  set1.intron.5.transition set1.intron.5.transition


##  Compute log-likelihood exon matrices, assuming intron
##  matrices describing background probabilities

gawk -f pro2log_ini.awk set1.intron.5.initial set1.cds.5.initial \
  >  set1.cds-intron.5.initial
gawk -f pro2log_tran.awk set1.intron.5.transition set1.cds.5.transition \
  >  set1.cds-intron.5.transition

gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$3}' \
  set1.cds-intron.5.initial > set1.cds-intron.5.initial.geneid
gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$4}' \
  set1.cds-intron.5.transition > set1.cds-intron.5.transition.geneid
