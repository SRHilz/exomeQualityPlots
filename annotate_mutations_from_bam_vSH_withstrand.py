#Notes: This script is an adaption of annotate_mutations_from_bam.py from the exome pipeline version 1.
#     It has been modified in the following ways:
#     1. The variable names 'rnaFilename','rnaString', and 'rna' have been changed to 'pileupFilename',
#             'pileupString', and 'pileup' to avoid confusion.
#     2. The pileup command has had the option '-s' added to its call (this adds an extra column of
#             output containing mapq scores)
#     3. The mapq scores, as well as quality scores, are now collected in getNormalCoverage() in the
#             arrays 'ref_qual','alt_qual','ref_mqual', and 'alt_mqual', and are included in the returns
#             of getNormalCoverage(). These are all output to the file qualFile, whose file handle is now passed
#             into annotate_mutations_from_bam(), and which is used to produce a graph for each varient in R.
#     4. A new variable, mqualityCutoff, has been added to the function getNormalCoverage(). This sets
#             cutoff for mapping quality.
#     5. To additionally filter for mapping quality in addition to sequenced base quality, the
#             statement 'if (ord(qualityString[qIndex])-33 >= qualityCutoff):' is now
#            'if ord(qualityString[qIndex])-33 >= qualityCutoff and ord(mqualityString[qIndex])-33 >= mqualityCutoff:'
#     6. To ensure that quality informaiton is not output more than once for every sample and variant, added the
#             array 'qualsAlreadyOutput', which, for every sample, keeps track of which loci the quality info has already
#             been output for.
#     7. To make more user-friendly, the variable 'bamPath' was added to get_samples_from_patient(),
#             which allows the user to specify the location of the Patient folders containing bamfiles are.
#             This path is also now an argument (sys.argv[5]) when calling annotate_mutations_from_bam_vSH.
#             The line 'col_bam = h.index("bam_file")', which was meant to obtain a path from the conversionfile
#             file, was also commented out. Finally, 'mutfile.split('.txt')[0]' was replaced with
#             'mutfile.split('/')[-1].split('.txt')[0], which makes it so that the output files go not to the directory
#             that contains the .snvs files, but rather, the current working directory.

import sys, subprocess, os

def get_samples_from_patient(mutfile, conversionfile, patient_ID, projectname, bamPath):
    ## read ID conversion file
    data = open(conversionfile).readlines()
    header = data[0]
    data = data[1:]
    ## parse header
    h = header.strip().split()
    col_pat = h.index("patient_ID")
    col_lib = h.index("lib_ID")
    col_st = h.index("sample_type")
    #col_bam = h.index("bam_file")

    ## pull out patient specific info
    data_pat = filter(lambda x:x.strip().split()[col_pat] == patient_ID, data)

    ## determine where files are stored
    #testID = data_pat[0].strip().split()[col_lib]
    #fullpath = "/data/jocostello/" + projectname + "/exomes/" + patient_ID + "/"
    #fileheader = ".bwa.realigned.rmDups"
    #fileheader = ".bwa.realigned.rmDups.recal"
    #if not os.path.isfile(fullpath + testID + fileheader + ".bam"):
    #    fullpath = fullpath.replace("data", "home", 1)
    #    if not os.path.isfile(fullpath + testID + fileheader + ".bam"):
    #        print "ERROR: files can not be found"
    #        print fullpath + testID + fileheader + ".bam"
    #        sys.exit(1)

    ## for each sample, call annotate_mutations_from_bam
    qualFile = open(patient_ID+'.qualityinfo.tmp','w')
    for line in data_pat:
        print mutfile
        sample = line.strip().split()[col_st]
        print sample

        #if patient_ID == "Patient126" and    sample != "Normal":
        #    fullpath = "/costellolab/jocostello/LG3/exomes/" + line.strip().split()[col_lib] + "/"
        #    fileheader = ".merged.sorted"

        #bamfile = fullpath + line.strip().split()[col_lib] + fileheader + ".bam"
        bamfile = bamPath+'/'+patient_ID+'/'+line.strip().split()[col_lib]+'.bwa.realigned.rmDups.recal.bam'
        print bamfile

        annotate_mutations_from_bam(mutfile, bamfile, sample, qualFile)
        mutfile = mutfile.split('/')[-1].split('.txt')[0] + '.%sQ.txt'%(sample)
        print mutfile
    qualFile.close()

    ## rename final mutfile
    command = ['mv', mutfile, patient_ID + ".snvs.anno.txt"]
    task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr)=task.communicate()


def annotate_mutations_from_bam(mutfile, bamfile, sample, qualFile):
    ## read mutation data
    data = open(mutfile).readlines()
    header = data[0]
    data = data[1:]

    ## parse header
    h = header.strip().split('\t')
    chr = h.index('contig')
    pos = h.index('position')
    ref = h.index('ref_allele')
    alt = h.index('alt_allele')
    alg = h.index('algorithm')

    ## update header
    h.append('%s_ref_reads'%(sample))#; r = h.index('%s_ref_reads'%(sample))
    h.append('%s_alt_reads'%(sample))#; a = h.index('%s_alt_reads'%(sample))
    h.append('%s_ref_Q20reads'%(sample))#; rQ = h.index('%s_Q20ref_reads'%(sample))
    h.append('%s_alt_Q20reads'%(sample))#; aQ = h.index('%s_Q20alt_reads'%(sample))

    ## prepare outfile
    annofile = open(mutfile.split('/')[-1].split('.txt')[0] + '.%sQ.txt'%(sample), 'w')
    annofile.write('\t'.join(h) + '\n')

    tmp_file_header = mutfile.split('/')[-1].split('.txt')[0] + '_' + sample

    ## generate bedfile for this patient
    print 'making bedfile'
    bedfile = tmp_file_header+ '.bed'
    outfile = open(bedfile, 'w')
    for m in data:
        split = m.split('\t')
        outfile.write(split[chr] + '\t' + split[pos] + '\t' + str(int(split[pos])+1) + '\n')
    outfile.close()

    ## generate tmp bam files split by read -f 0x80 Z00371.bwa.realigned.rmDups.recal.bam | wc
    bamfileR1 = tmp_file_header + '.r1.bam'
    command = ['/home/jocostello/tools/samtools-0.1.12a/samtools', 'view', '-f','0x40', '-b', '-o', bamfileR1, bamfile]
    task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr)=task.communicate()
    print stderr
    bamfileR2 = tmp_file_header + '.r2.bam'
    command = ['/home/jocostello/tools/samtools-0.1.12a/samtools', 'view', '-f','0x80', '-b', '-o', bamfileR2, bamfile]
    task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout,stderr)=task.communicate()
    print stderr

    ## generate Read1 and Read2 mpileup for this patient
    print 'making pileup'
    mpilefile = tmp_file_header + '.pileup'
    outfile = open(mpilefile, 'w')
    dic = {bamfileR1 : '1', bamfileR2 : '2'}
    pileupdic = {}
    for splitbam in [bamfileR1,bamfileR2]:
        command = ['/home/jocostello/tools/samtools-0.1.12a/samtools', 'pileup', '-s','-l', bedfile, splitbam]
        task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        (stdout,stderr)=task.communicate()
        pileupout = stdout.split('\n')
        for line in pileupout:
            line = line.rstrip().split('\t')
            if len(line) == 7:
                if line[0] not in pileupdic:
                    pileupdic[line[0]] = {}
                if line[1] not in pileupdic[line[0]]:
                    pileupdic[line[0]][line[1]] = line
                    pileupdic[line[0]][line[1]].append(dic[splitbam]*int(line[3]))
                else:
                    pileupdic[line[0]][line[1]][3] = int(pileupdic[line[0]][line[1]][3]) + int(line[3])
                    pileupdic[line[0]][line[1]][4] = pileupdic[line[0]][line[1]][4] + line[4]
                    pileupdic[line[0]][line[1]][5] = pileupdic[line[0]][line[1]][5] + line[5]
                    pileupdic[line[0]][line[1]][6] = pileupdic[line[0]][line[1]][6] + line[6]
                    pileupdic[line[0]][line[1]][7] = pileupdic[line[0]][line[1]][7] + dic[splitbam]*int(line[3])
            elif len(line) >1:
                print line
                print 'Pileup output does not have the expected format. Please check pileup version. Press ^c to exit.'
    for chromosome in pileupdic:
        print pileupdic[chromosome].keys()
        for position in pileupdic[chromosome]:
                outfile.write('\t'.join([str(x) for x in pileupdic[chromosome][position]])+'\n')
    outfile.close()

    print 'annotating coverage'
    qualsAlreadyOutput = []#as there are multiple entries for the same positin in the mutfile and thus bed, this ensures that qual information is only output once per locus per sample
    for line in data:
        l = line.strip('\n').split('\t')
        if l[alg] != "MuTect": (countR, countA, countQR, countQA) = ("NA", "NA", "NA", "NA")
        else:
            (countR, countA, countQR, countQA, qualR, qualA, mqualR, mqualA, strandR, strandA, readR, readA) = getNormalCoverage(mpilefile, l[chr], l[pos], l[ref], l[alt], 20, 0)
            variantID = l[chr]+'_'+l[ref]+l[pos]+l[alt]
            if variantID not in qualsAlreadyOutput:
                for x in range(0,len(qualR)):
                    qualFile.write(str(qualR[x])+'\t'+str(mqualR[x])+'\tref\t'+sample+'\t'+variantID+'\t'+strandR[x]+'\t'+readR[x]+'\n')
                for x in range(0,len(qualA)):
                    qualFile.write(str(qualA[x])+'\t'+str(mqualA[x])+'\talt\t'+sample+'\t'+variantID+'\t'+strandA[x]+'\t'+readA[x]+'\n')
                qualsAlreadyOutput.append(variantID)
        l.append(str(countR))
        l.append(str(countA))
        l.append(str(countQR))
        l.append(str(countQA))
        annofile.write('\t'.join(l) + '\n')

    ## clean up
    print 'cleaning up'
    for f in [bedfile, bamfileR1, bamfileR2, mpilefile]:
        command = ['rm' , '-f', f]
        task=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        (stdout,stderr)=task.communicate()

    annofile.close()


def getNormalCoverage(pileupFilename, contig, position, ref_allele, alt_allele, qualityCutoff = 0, mqualityCutoff = 0):
    """Parses a mpileup file and returns the normal coverage stats"""
    ref_count = 0
    alt_count = 0
    ref_mQcount = 0
    alt_mQcount = 0
    pileupString = ""
    ref_qual = []
    alt_qual = []
    ref_mqual = []
    alt_mqual = []
    ref_strand = []
    alt_strand = []
    ref_read = []
    alt_read = []
    skip = 0

    pileup = open(pileupFilename)

    for line in pileup:
        line = line.rstrip("\n")
        line = line.split("\t")
        if (line[0] == contig):
            if (line[1] == position):
                pileupString = line[4]
                qualityString = line[5]
                mqualityString = line[6]
                readString = line[7]
                break
    pileup.close()

    #print pileupString, qualityString

    #process the reads from the pileup
    qIndex = 0
    for p in xrange(0,len(pileupString)):
        base = pileupString[p]
        if (skip > 0):
            skip -= 1
            continue
        if (base == "^"):    # this marks the beginning of a read, the next position is the read mapping quality
            skip = 1
            continue
        if (base == "+" or base == "-"):    ## indicates that what follows is an indel
            continue
        if (base.isdigit() == True):    ## indicates length of an indel, the following positions are the indel sequence
            if (pileupString[p+1].isdigit()): skip = int(pileupString[p:p+2])+1
            else: skip = int(base)
            #print skip
            continue
        if (base == "$"):    ## this marks the end of the read
            continue
        ## none of the above have an associated value in the quality string, so only increment qIndex after this
        if base.upper() == ref_allele:
            ref_count += 1
            ref_qual.append(ord(qualityString[qIndex])-33)
            ref_mqual.append(ord(mqualityString[qIndex])-33)
            ref_read.append(readString[qIndex])
            if ord(qualityString[qIndex])-33 >= qualityCutoff and ord(mqualityString[qIndex])-33 >= mqualityCutoff:
                ref_mQcount +=1
            if base.isupper():
                ref_strand.append('+')
            else: ref_strand.append('-')

        elif base.upper() == alt_allele:
            alt_count += 1
            alt_qual.append(ord(qualityString[qIndex])-33)
            alt_mqual.append(ord(mqualityString[qIndex])-33)
            alt_read.append(readString[qIndex])
            if ord(qualityString[qIndex])-33 >= qualityCutoff and ord(mqualityString[qIndex])-33 >= mqualityCutoff:
                alt_mQcount +=1
            if base.isupper():
                alt_strand.append('+')
            else: alt_strand.append('-')
        ## increment qIndex
        qIndex +=1

    return (ref_count, alt_count, ref_mQcount, alt_mQcount, ref_qual, alt_qual, ref_mqual, alt_mqual, ref_strand, alt_strand, ref_read, alt_read)



if __name__=="__main__":
    if len(sys.argv) != 6:
        print 'usage: %s mutationfile.txt patient_ID_conversions.txt patient_ID projectname bamPath' %(sys.argv[0])
        sys.exit(1)

    get_samples_from_patient(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].strip(), sys.argv[5].strip())
