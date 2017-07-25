#!/bin/bash
cd $PBS_O_WORKDIR
bampath="/costellolab/jocostello/LG3/exomes_recal/${patient}/" 
#bampath="/home/shilz/data/exome2_recal/${patient}/"
for entry in $bampath*.bam
do
    output=$(echo "$entry" | sed "s|$bampath||g").hist
    bedtools coverage -hist -b $entry -a /home/shilz/database/nimblegen3/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed | grep ^all > $output 
    done
