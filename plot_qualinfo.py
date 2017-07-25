#Stephanie Hilz
#Costello Lab
#2016.10.27

#Usage: python plot_qualinfo.py <hg19.faFilePath> <R.mutationsFilePath> <qualinfo.tmpFilePath>
#Desription: Cross-references the qualinfo.tmp file, which has quality info for each variant detected
#  by Mutect, with the R.mutations file for that patient, which only has the final filtered set of
#  variants. The output is a final qualinfo.txt file that is then read into R to produce the final
#  qual figure for each variant.
#Input:
#  hg19.faFilePath - path to fasta for ref genome used
#  R.mutationsFilePath - this is a tab delim txt file which lists the final, filtered mutations
#    for a given patient.
#  qualinfo.tmpFilePath - this is a tab delim file which lists quality information for each sample
#    at each site where a variant was detected

#Dependencies
import sys
from Bio import SeqIO

#Funcitons
def getcontext(SeqIOdic,chromosome, position, upstream, downstream):
    #Usage: Retrieves context information of a specified size for centered at chromosome and position.
    #Input:
        # <SeqIOdic> the output of SeqIO.to_dict(SeqIO.parse(open(faFile),'fasta')) for a given faFile (ex. hg19.fa)
        # <chromosome> column in coordFile containing chromosome information (should correspond to .fa file)
        # <position> column in coordFile containing position information (should correspond to .fa file)
        # <upstream> number of bases you want upstream of the genomic coordinate
        # <downstream> number of bases you want downstream of the genomic coordinates

    if position.isdigit():
        position = int(position)-1#essential to subtract one so index is base 0
        if chromosome in SeqIOdic:
            if position < len(SeqIOdic[chromosome]):
                contextChunk = ''.join(SeqIOdic[chromosome][position-int(upstream):position+int(downstream)+1])
            else: print 'Position specified out of range for chromosome: recheck chromosome and position'
        else: print 'Chromosome specified is not in .fa file'
    else: print 'Position specified is not an integer. Please enter an integer.'

    return contextChunk.upper()

#User-specified variables:
qualCutoff = 20
mqualCutoff = 0
calledAgainst = 'Normal'

#Part 1: Import and index fasta
print 'Part 1: Reading in ',sys.argv[1]
faFile = sys.argv[1]
record_dic = SeqIO.to_dict(SeqIO.parse(open(faFile),'fasta'))

#Part 2: Get list of final mutants for patient + lookup trinucleotide content + output var info
## open R.mutations file
print 'Part 2: Getting final called variant information from ',sys.argv[2]
mutationFile = sys.argv[2]
data = open(mutationFile, 'rU').readlines()
header = data[0]
data = data[1:]
## parse header
h = header.strip().split('\t')
gene = h.index('#gene')
chr = h.index('contig')
pos = h.index('position')
ref = h.index('ref_allele')
alt = h.index('alt_allele')
alg = h.index('algorithm')
## pull out mutation list
dic = {}#will hold all MuTect variants, indexed by contig_refpositionalt ID, giving gene_contig_refpositionaltID
vardic = {}#will hold all MuTect variant info, indexed by IDplus
dtoxogdic = {}#will hold all Mutext variant info needed for D-ToxoG, indexed by IDplus
for line in data:
    l = line.rstrip().split('\t')
    if l[alg] == "MuTect":
        ID = l[chr]+'_'+l[ref]+l[pos]+l[alt]#unique base specific ID
        IDplus = l[gene]+'_'+ID
        dic[ID] = IDplus
        tricontext = getcontext(record_dic,l[chr], l[pos], 1, 1)
        dtoxogdic[IDplus] =[l[chr],l[pos],l[pos],l[ref],l[alt],'','',tricontext,0,0,0,0,0,'SNP']
        tricontext_collapsed = tricontext[0]+'N'+tricontext[2]
        vardic[IDplus] = [l[ref],l[alt],tricontext_collapsed,l[-1]]
        outFile = open('varinfo.txt','w')
for IDplus in vardic:
    outFile.write(IDplus+'\t'+'\t'.join(vardic[IDplus])+'\n')
outFile.close()
print 'Variant info output to varinfo.txt'

#Part 3: Filter qualinfo file and output
## prepare outfile
print 'Part 3: Getting quality information from ',sys.argv[3]
dtoxog_counts = {}
qualinfoFile = sys.argv[3]
patientID = qualinfoFile.split('/')[-1].split('.')[0]
outFile = open(patientID+'.qualityinfo.txt','w')
## filter qualinfo.tmp file and output
inFile = open(qualinfoFile,'rU')
while 1:
    line = inFile.readline()
    if not line: break
    sline = line.rstrip().split('\t')
    if sline[4] in dic:
        sline[4] = dic[sline[4]]
        sample = sline[3]
        variant = sline[4]
        allele = sline[2]
        read = str(sline[6])
        strand = sline[5]
        called = 'N'
        calls = vardic[variant][3].split(',')
        for item in calls:
            if item == sample:
                called = 'Y'
        outFile.write('\t'.join(sline)+'\t'+'\t'.join(vardic[variant][:3])+'\t'+called+'\n')
        if int(sline[0]) >= qualCutoff and int(sline[1]) >= mqualCutoff:
            if sample not in dtoxog_counts:
                dtoxog_counts[sample] = {}
            if variant not in dtoxog_counts[sample]:
                dtoxog_counts[sample][variant] = {}
            if allele not in dtoxog_counts[sample][variant]:
                dtoxog_counts[sample][variant][allele] = {}
            if read not in dtoxog_counts[sample][variant][allele]:
                dtoxog_counts[sample][variant][allele][read]={}
            if strand not in dtoxog_counts[sample][variant][allele][read]:
                dtoxog_counts[sample][variant][allele][read][strand]=0
            dtoxog_counts[sample][variant][allele][read][strand]+=1 #will be counts of [sample][var][alt/ref][r1/r2]
inFile.close()
outFile.close()
print 'Final quality info output to qualityinfo.txt'

# #Part 4: Create D-ToxoG input file
# print 'Part 4: Creating D-ToxoG output file'
# outFile = open('dtoxog_input.txt','w')
# outFile.write('Chromosome\tStart_position\tEnd_position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tref_context\ti_t_ALT_F1R2\ti_t_ALT_F2R1\ti_t_REF_F1R2\ti_t_REF_F2R1\ti_t_Foxog\tVariant_Type\n')
# for sample in dtoxog_counts:
#     for variant in dtoxog_counts[sample]:
#         ref_F1R2 = 0
#         ref_F2R1 = 0
#         alt_F1R2 = 0
#         alt_F2R1 = 0
#         local_dtoxoglist = dtoxogdic[variant][:]#creates a local copy of list for variant that we can modify without changing dtoxogdic
#         for allele in dtoxog_counts[sample][variant]:
#             for read in dtoxog_counts[sample][variant][allele]:
#                 for strand in dtoxog_counts[sample][variant][allele][read]:
#                     if allele == 'ref':
#                         if (read == '1' and strand == '+') or (read == '2' and strand == '-'):
#                             ref_F1R2 += dtoxog_counts[sample][variant][allele][read][strand]
#                         else:
#                             ref_F2R1 += dtoxog_counts[sample][variant][allele][read][strand]
#                     if allele == 'alt':
#                          if (read == '1' and strand == '+') or (read == '2' and strand == '-'):
#                              alt_F1R2 += dtoxog_counts[sample][variant][allele][read][strand]
#                          else:
#                              alt_F2R1 += dtoxog_counts[sample][variant][allele][read][strand]
#         local_dtoxoglist[5] = sample
#         local_dtoxoglist[6] = calledAgainst
#         local_dtoxoglist[8:12] = [alt_F1R2, alt_F2R1, ref_F1R2, ref_F2R1]
#         if alt_F1R2 > 0 or alt_F2R1 > 0:
#             if local_dtoxoglist[3] == 'C' or local_dtoxoglist[3] == 'A':
#                 foxog = float(alt_F2R1)/(alt_F1R2 + alt_F2R1)
#             else:
#                 foxog = float(alt_F1R2)/(alt_F1R2 + alt_F2R1)
#             local_dtoxoglist[12] = foxog
#             outFile.write('\t'.join([str(x) for x in local_dtoxoglist])+'\n')
# outFile.close()
