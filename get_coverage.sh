#!/bin/bash                                                                                                                                                                                                        
cd $PBS_O_WORKDIR
bampath="/costellolab/jocostello/LG3/exomes_recal/${patient}/"
#bampath="/home/shilz/data/exome2_recal/${patient}/"                                                                                                                                                               
for entry in $bampath*.bam
do
    output=$(echo "$entry" | sed "s|$bampath||g").hist
    bedtools coverage -hist -b $entry -a $bedpath | grep ^all > $output
    done
