#!/bin/bash
#
##
#$ -clear
#$ -S /bin/bash
#$ -cwd
#$ -j y
#

cd $PBS_O_WORKDIR

## run annotation code
python ${scriptpath}/annotate_mutations_from_bam_vSH_withstrand.py ${mutationfile} ${conversionfile} ${patient} ${project} ${bampath}

## remove intermediate files
rm -f ${patient}.snvs.*Q*.txt

echo "Finished"
