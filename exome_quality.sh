#!/bin/bash
#
##
#$ -clear
#$ -S /bin/bash
#$ -cwd
#$ -j y
#

cd $PBS_O_WORKDIR

## reformat quality info file
python ${scriptpath}/plot_qualinfo.py ${reffasta} ${mutationfile} ${patient}.qualityinfo.tmp

## make plots
Rscript --vanilla ${scriptpath}/plot_qualinfo.R ${patient} ${patient}.qualityinfo.txt ${mutationfile}

## make coverage histograms
Rscript --vanilla ${scriptpath}/coveragePlots.R ${conversionfile} ${patient} ./ ${patient}_qualplots/libraryQuality

## make TERT coverage plots
python ${scriptpath}/afTERThought.py ${bampath}/${patient} ${conversionfile} ${patient} ${patient}_qualplots/VAFPatterns

## remove intermediate files
#rm -f ${Patient}.qualityinfo.tmp

echo "Finished"
