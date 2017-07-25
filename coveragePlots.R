# Name: coveragePlots
# Creator: Stephanie Hilz (UCSF Costello Lab)
# Date: 2017.01.12
# Description: Creates plots to visualize exome coverage based on bedtools coverage
#  histogram files. Adapted from:http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
# Usage: R coveragePlots <patient_ID_conversions.txt> <patient ID> <bedtools file path> <outputdir>
# Input: 
#  <patient_ID_conversion.txt> path for Tali's sample ID conversion file
#  <patient ID> ID for patient (ex. 'Patient300')
#  <bedtools file path> path to directory containing file(s) produced by bedtools (command to produce these files is something like bedtools coverage -hist -b $entry -a /home/shilz/database/nimblegen3/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed | grep ^all > $output)
#  <output dir> name of directory to output to
# Output: 
#  <PatientID.coverage.png> plot of exome coverage distribution for all samples for that patient

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Usage: R coveragePlots <patient_ID_conversions.txt> <patient ID> <bedtools file path> <output dir>", call.=FALSE)
}

# Set variables
conversionfile <- args[1]
patientID <- args[2]
filepath <- args[3]
outputdir <- args[4]

# Get a list of the bedtools output files you'd like to read in
print(filesUnordered <- list.files(path=filepath, pattern=".hist"))

# Reorder files (requires patient_ID_conversions.txt file)
print('here')
minLast <- min(count.fields(conversionfile, sep = "\t")) #added/edited next few lins 2017.06.21 after someone messed with conversion file
maxLast <- max(count.fields(conversionfile, sep = "\t"))
conversion <- read.csv(conversionfile,header=FALSE,colClasses=c(rep("character", minLast), rep("NULL",maxLast-minLast)), sep='\t',stringsAsFactors = FALSE)
colnames(conversion) <- conversion[1,]#just leaving in original col name as won't cause probs, grep out what I want from file
labs <- sort(conversion[which(conversion$patient_ID==patientID),]$sample_type)
files <- c()
for (element in labs){
  element <- gsub('\\.','-',element)
  sampleID <- conversion[which(conversion$patient_ID==patientID & conversion$sample_type==element),]$lib_ID
  files <- append(files, filesUnordered[grep(sampleID, filesUnordered)])
}

print('here')
# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(paste(filepath,'/',files[i],sep=''))
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

# Pick some colors
cols <- c("#A6CEE3","#1F78B4",'#696969',"#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",'#00F5FF','#556B2F','#7CFC00','#8B0A50','#000080','#FFC125','#D3D3D3')

# Save the graph to a file
pdf(paste(outputdir,"/coverage.pdf",sep=''))

# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
legend("topright", legend=labs, col=cols, lty=1, lwd=4)

dev.off()
