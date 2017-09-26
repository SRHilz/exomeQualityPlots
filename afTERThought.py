#Stephanie Hilz
#Costello Lab
#2017.01.16

#Usage: python afTERThought.py <path to directory with bams> <path to conversionfile> <Patient ID> <outputdir>
#Desription: Pulls out any reads corresponding to the 1295228 and 1295250 TERT promoter mutations
#   and outputs a graph showing the number of reads supporting WT, G>A, and other
#Input:
#  <path to directory with bams> - this is the path to the directory containing the bams. Will
#    all bams in that directory
#  <path to conversionfile> - this is the path to Tali's conversion file. This is just so I can
#    use logical names for the sample names
#  <Patient ID> - ex: Patient300
#  <outputdir> - name of output directory


#Dependencies
import pysam
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import fnmatch

#Part 1: Reading in conversion file
print 'Part 1: Reading in conversion file.'
conversionfile = sys.argv[2]
patientID = sys.argv[3]
conversion = open(conversionfile,'r')
conversionDic = {}
while 1:
    line = conversion.readline()
    if not line: break
    sline = line.rstrip().split('\t')
    if sline[2] == patientID:
        conversionDic[sline[3]] = sline[0]
conversion.close()
orderedNames = sorted(conversionDic.keys())
print orderedNames

#Part 2: Count G, A, and other reads
print 'Part 2: Reading in .bam files.'
counts = {}
chrom = 'chr5'
loci = ['1295228','1295250']
bases = ['G','A','C or T']
for locus in loci:
    counts[locus] = {}
    for base in bases:
        counts[locus][base] = {}
dirPath = sys.argv[1]
i = 0 #will be index for files
for file in os.listdir(dirPath):
    if fnmatch.fnmatch(file, '*.bam.bai'):
        print file
        bamfile = file.replace('.bai','')
        samfile = pysam.Samfile(dirPath+'/'+bamfile,'rb')
        sampleID = bamfile.split('.')[0]
        if 'trim' in sampleID:#added 2017.09.26; specific for Batch16
            sampleID = sampleID.replace('-trim','')
        for locus in loci:
            for base in bases:
                counts[locus][base][sampleID] = 0
            for pileupcolumn in samfile.pileup(chrom, int(locus)-1, int(locus)):
                for pileupread in pileupcolumn.pileups:
                    if pileupcolumn.pos == int(locus)-1:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == 'C' or base == 'T': base = 'C or T'
                        counts[locus][base][sampleID]+=1
print counts
#Part 3: Plot
print 'Part 3: Plotting.'
outputdir = sys.argv[4]
index = np.arange(len(orderedNames))  # the x locations for the groups
width = 0.2       # the width of the bars
color = ['r','b','g']
plt.figure(1)
plt.gcf().subplots_adjust(bottom=0.2)
i = 0
while i <len(loci):
    maximumValue = 0
    plt.subplot(211+i)
    rects = []
    j = 0
    while j<len(bases):
        toPlot = []
        k = 0
        while k<len(orderedNames):
            value = counts[loci[i]][bases[j]][conversionDic[orderedNames[k]]]
            toPlot.append(value)
            if value > maximumValue: maximumValue = value
            k+=1
        rects.append(plt.bar(index+(j*width), toPlot, width, color=color[j], label=bases[j]))
        j+=1
    plt.title(loci[i])
    plt.legend()
    plt.ylabel('Counts')
    plt.ylim((0,maximumValue+5))
    if (i==1):
        plt.xticks(index + width, orderedNames, rotation='vertical')
    else:
        plt.xticks(index + width, ['']*len(orderedNames), rotation='vertical')
    i+=1

plt.savefig(outputdir+'/TERTp.pdf')
