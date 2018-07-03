# exomeQualityPlots

The current setup is to run code in three steps per patient.

#### 1. Generate coverage distribution files, which are used in Step #3 to generate read coverage plots for each library.

The script used for this step is:

```
get_coverage.sh
```

Example of usage with the required arguments (note: can require alot of memory, so if fails, rerun with more memory):

```
patient=Patient300

qsub -l vmem=300gb -v patient=$patient,bedpath=/home/shilz/database/exomeKitBeds/SeqCap_EZ_Exome_v3_capture.bed /home/shilz/tools/exomeQualityPlots/get_coverage.sh
```

If successful, should output a file for each library (including normal) of non-zero size with the extension .hist.

#### 2. Collect quality information for each position where a mutation was called, which are used in Step #3 to generate quality plots. This step can be performed in parallel with Step 1.

The script used for this step is:

```
annotate_mutations_from_bam_vSH.sh
```

Example of usage with the required arguments:

```
patient=Patient300

qsub -l vmem=96gb -v mutationfile=/costellolab/jocostello/LG3/MutInDel/$patient.snvs.indels.filtered.overlaps.txt,conversionfile=/costellolab/mazort/LG3/exome/patient_ID_conversions.txt,patient=$patient,project=LG3,bampath=/costellolab/jocostello/LG3/exomes_recal/,scriptpath=/home/shilz/tools/exomeQualityPlots/ /home/shilz/tools/exomeQualityPlots/annotate_mutations_from_bam_vSH.sh
```

If successful, should output a file called ${patient}.qualityinfo.tmp

#### 3. Make final quality plots using files generated in Steps 1 and 2.

The script used for this step is:

```
exome_quality.sh
```

Example of usage with the required arguments:

```
patient=Patient300

qsub -l vmem=96gb -v scriptpath=/home/shilz/tools/exomeQualityPlots/,reffasta=/home/jocostello/shared/LG3_Pipeline/resources/UCSC_HG19_Feb_2009/hg19.fa,patient=$patient,mutationfile=/costellolab/jocostello/LG3/MutInDel/$patient.snvs.indels.filtered.overlaps.txt,conversionfile=/costellolab/mazort/LG3/exome/patient_ID_conversions.txt,bampath=/costellolab/jocostello/LG3/exomes_recal/ ~/tools/exomeQualityPlots/exome_quality.sh
```

If successful, should output the following:
- ${patient}.qualityinfo.txt
- ${patient}.qualitystats.txt
- a directory called ${pateint}_qualplots/, with the following subdirectories containing the following files:
  - variantQuality/
    - .pdfs for each mutation called by MuTect for all patients with less than <300 mutations called
  - variantSpectra/ (will not output all of these if >300 mutations called for the patient)
    - allReads_counts.pdf
    - allReads_proportion.pdf
    - allVariantsPresent_counts.pdf
    - allVariantsPresent_proportion.pdf
    - allVariantsPresent_trinucleotide_counts.pdf
    - allVariantsPresent_trinucleotide_proportion.pdf
    - calledVariants_counts.pdf
    - calledVariants_proportion.pdf
    - typefrequency.pdf
  - libraryQuality/
    - coverage.pdf
    - overall_basequality.pdf
    - overall_mappingquality.pdf
    - ref_v_alt_basequality.pdf
    - rev_v_alt_mappingquality.pdf
  - VAFPatterns/
    - calldistribution.pdf
    - commondrivers.pdf
    - pca.pdf
    - TERTp.pdf
    
