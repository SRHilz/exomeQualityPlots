#dependencies
if("dplyr" %in% rownames(installed.packages()) == FALSE) {
  install.packages("dplyr", repos="http://cran.rstudio.com/")
  }
if("plyr" %in% rownames(installed.packages()) == FALSE) {
  install.packages("plyr", repos="http://cran.rstudio.com/")
}
if("scales" %in% rownames(installed.packages()) == FALSE) {
  install.packages("scales", repos="http://cran.rstudio.com/")
}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2", repos="http://cran.rstudio.com/")
}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {
  install.packages("RColorBrewer", repos="http://cran.rstudio.com/")
}
if("gtools" %in% rownames(installed.packages()) == FALSE) {
  install.packages("gtools", repos="http://cran.rstudio.com/")
}
library(dplyr)
library(plyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(gtools)

extractMuts <- function(variantID){
  s <- strsplit(variantID, "_")
  s <- unlist(s)
  s <- s[3]
  ref_base <- substr(s,1,1)
  alt_base <- substr(s,nchar(s),nchar(s))
  mut <- paste(ref_base,'>',alt_base, sep='')
  return(mut)
}

getRatio <- function(x){
  toReturn <- x/sum(x)
  return(toReturn)
}

getVAFs <- function(input, samples, qualCutoff, mqualCutoff){
  detach("package:dplyr", unload=TRUE)
  qual_cutoff <- qualCutoff
  mqual_cutoff <- mqualCutoff
  data_filtered <- input[which(input$V1>=qual_cutoff & input$V2>=mqual_cutoff),]
  data_filtered$V1 <- NULL
  data_filtered$V2 <- NULL
  data_filtered$V6 <- NULL
  colnames(data_filtered) <- c('type','sample','variant')
  variantIDs <- unique(data_filtered$variant)
  data_filtered <- count(data_filtered, vars = c("type","sample","variant"))
  VAFs_df <- c()
  for (i in 1:length(variantIDs)){
    variantVAFs <- c()
    for (j in 1:length(samples)){
      refCount <- 0
      altCount <- 0
      if (length(data_filtered[which(data_filtered$type == 'ref' & data_filtered$sample == samples[j] & data_filtered$variant == variantIDs[i]),]$freq) > 0){
        refCount <- data_filtered[which(data_filtered$type == 'ref' & data_filtered$sample == samples[j] & data_filtered$variant == variantIDs[i]),]$freq
      }
      if (length(data_filtered[which(data_filtered$type == 'alt' & data_filtered$sample == samples[j] & data_filtered$variant == variantIDs[i]),]$freq) > 0){
        altCount <- data_filtered[which(data_filtered$type == 'alt' & data_filtered$sample == samples[j] & data_filtered$variant == variantIDs[i]),]$freq
      }
      variantVAFs <- append(variantVAFs, altCount/(altCount + refCount))
    }
    VAFs_df <- rbind(VAFs_df, variantVAFs)# adds a row containing VAFs for all samples for that varient
  }
  colnames(VAFs_df) <- samples
  rownames(VAFs_df) <- variantIDs
  return (VAFs_df)
}

# OVERALL SETUP and VARIABLE PARAMETERS
# Args supplied at command line
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Usage: Rscript --vanilla plot_qualinfo.R <patient ID> <path to quality info file> <path to mutations.R file>", call.=FALSE)
}
patientID <-args[1] 
qualinfofile <-args[2]# paste0('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Variants/',patientID,'.qualityinfo.txt')
 # #
mutationsfile <- args[3]#paste0('/Users/srhilz/Documents/Professional/Positions/UCSF_Costello/Data/',patientID,'/Variants/',patientID,'.R.mutations.txt') # 

# Hardcoded settings ordering mutational spectra output
substitutions <- c('A>C','T>G','A>G','T>C','A>T','T>A','C>G','G>C','C>A','G>T','C>T','G>A')
substitutionsColor <- brewer.pal(12,'Set3')
substitutionsColorset <- c('A>C'=substitutionsColor[1],'T>G'=substitutionsColor[2],'A>G'=substitutionsColor[3],'T>C'=substitutionsColor[4],'A>T'=substitutionsColor[5],'T>A'=substitutionsColor[6],'C>G'=substitutionsColor[7],'G>C'=substitutionsColor[8],'C>A'=substitutionsColor[9],'G>T'=substitutionsColor[10],'C>T'=substitutionsColor[11],'G>A'=substitutionsColor[12])
substitutions_grouping <- c('A:T>C:G','A:T>C:G','A:T>G:C','A:T>G:C','A:T>T:A','A:T>T:A','C:G>G:C','C:G>G:C','C:G>A:T','C:G>A:T','C:G>T:A','C:G>T:A')
substitutionsGroupedColor <- brewer.pal(6,'Set3')
substitutionsGroupedColorset <- c('A:T>C:G'=substitutionsGroupedColor[1],'A:T>G:C'=substitutionsGroupedColor[2],'A:T>T:A'=substitutionsGroupedColor[3],'C:G>G:C'=substitutionsGroupedColor[4],'C:G>A:T'=substitutionsGroupedColor[5],'C:G>T:A'=substitutionsGroupedColor[6])
trinucleotide_grouping <- c('ANA','ANC','ANG','ANT','CNA','CNC','CNG','CNT','GNA','GNC','GNG','GNT','TNA','TNC','TNG','TNT')
trinucleotideColor <- colorRampPalette(brewer.pal(12, "Set3"))(16)
trinucleotideColorset <- c('ANA'=trinucleotideColor[1],'ANC'=trinucleotideColor[2],'ANG'=trinucleotideColor[3],'ANT'=trinucleotideColor[4],'CNA'=trinucleotideColor[5],'CNC'=trinucleotideColor[6],'CNG'=trinucleotideColor[7],'CNT'=trinucleotideColor[8],'GNA'=trinucleotideColor[9],'GNC'=trinucleotideColor[10],'GNG'=trinucleotideColor[11],'GNT'=trinucleotideColor[12],'TNA'=trinucleotideColor[13],'TNC'=trinucleotideColor[14],'TNG'=trinucleotideColor[15],'TNT'=trinucleotideColor[16])
qual_cutoff <- 20
mqual_cutoff <- 0

# Create master output directory
dir.create(paste(patientID,'_qualplots',sep=''),showWarnings = FALSE)

# Read in data and produce variant list and samples list
data <- read.table(qualinfofile, header=FALSE)  
samples <- mixedsort(as.character(unique(data$V4)))
variants <- as.vector(unique(data$V5))
data$V4 <- factor(data$V4, levels = as.factor(samples))#orders samples by mixedsort

# Read in mutations.R file, that will be used to generate the annotated mutations.R file
mutationdata <- read.table(mutationsfile, header=TRUE, sep='\t', comment.char="")
mutationdata <- mutationdata[which(mutationdata$algorithm=='MuTect'),]
mutationdata <- mutationdata[,c(1:9,dim(mutationdata)[2])]
rownames(mutationdata) <- paste(mutationdata$X.gene,'_',mutationdata$contig,'_',mutationdata$ref_allele,mutationdata$position,mutationdata$alt_allele,sep='')
mutationdata$fewalt <- 0 #a value of 1 denotes that there are less than 4 high-quality reads supporting the alt allele in all libraries
mutationdata$SB <- NA #a value of 0 indicates that for libraries with >=4 reads, there was at least one library where <=2x reads occured (no SB) from one strand vs the other; a value of 1 means there were no libraries supporting the varient without SB
mutationdata$OB <- NA #a value of 0 indicates that for libraries with >=4 reads, there was at least one library where <=2x reads occured (no OB) from one pair orientation vs the other; a value of 1 means there were no libraries supporting the varient without OB
mutationdata$qflag <- NA #a value of 0 indicates that there was at least one library where the base quality score distribution for reads supporting the variant allele was not different than those supporting the reference; a value of 1 means there were no libraries supporting the varient without a difference
mutationdata$mqflag <- NA #a value of 0 indicates that there was at least one library where the mapping quality score distribution for reads supporting the variant allele was not different than those supporting the reference; a value of 1 means there were no libraries supporting the varient without a difference

#VARIANT QUALITY DISTRIBUTION PLOTS
# create output dir
if (length(variants) <= 300){#this makes sure this section of code, which is computationally costly for lots of variants, is not executed for hypermutators
  dir.create(paste(patientID,'_qualplots/variantQuality',sep=''), showWarnings <- FALSE)
  data_filtered <- data[which(data$V1>=qual_cutoff & data$V2>=mqual_cutoff),]
  data_collapsed_strand <- count(data_filtered[,3:6], vars = c("V3","V4","V5","V6"))
  data_collapsed_orientation <- count(data_filtered[,3:7], vars = c("V3","V4","V5","V6","V7"))
  # create additional collapsed forms of data (counts of alt and ref for each vairant + sample)
  data_allcounts <- count(data[,3:5], vars = c("V3","V4","V5"))
  data_highqualcounts <- count(data_filtered[,3:5], vars = c("V3","V4","V5"))
  
  # set plot parameters
  plot_nrows = 4
  if(length(samples)%%4 == 0){
    plot_ncols = length(samples)/4
  } else {plot_ncols = floor(length(samples)/4)+1}
  #plot
  for (i in 1:length(variants)){
    allaltcounts <- 0
    pdf(paste(patientID,"_qualplots/variantQuality/",variants[i],".pdf", sep=""))
    par(mfrow=c(plot_nrows, plot_ncols),mar=rep(2,4))
    for (j in 1:length(samples)){
      plot(1, type="n", xlab="", ylab="",  xlim=c(0,70), ylim=c(0, .5), main=samples[j])
      if (length(data[which(data$V3=='ref' & data$V4==samples[j] & data$V5==variants[i]),]$V1)>1){
        refqual <- data[which(data$V3=='ref' & data$V4==samples[j] & data$V5==variants[i]),]$V1
        refmqual <- data[which(data$V3=='ref' & data$V4==samples[j] & data$V5==variants[i]),]$V2
        lines(density(refqual), col='firebrick', lty=1)
        lines(density(refmqual), col='blue',lty=1)
      }
      if (length(data[which(data$V3=='alt' & data$V4==samples[j] & data$V5==variants[i]),]$V1)>1){
        altqual <- data[which(data$V3=='alt' & data$V4==samples[j] & data$V5==variants[i]),]$V1
        altmqual <- data[which(data$V3=='alt' & data$V4==samples[j] & data$V5==variants[i]),]$V2
        lines(density(altqual), col='plum1',lty=1)
        lines(density(altmqual), col='deepskyblue',lty=1)
      }
      if (length(data[which(data$V3=='ref' & data$V4==samples[j] & data$V5==variants[i]),]$V1)>1 & length(data[which(data$V3=='alt' & data$V4==samples[j] & data$V5==variants[i]),]$V1)>1){
        p <- wilcox.test(refqual,altqual)$p.value
        qualpvalue <- paste('qual p=', round(p,3), sep='')
        if(is.na(mutationdata[variants[i],]$qflag) & is.finite(p)){
          if(p<0.05){#flag
            mutationdata[variants[i],]$qflag = 1
          } else {#mark as okay
            mutationdata[variants[i],]$qflag = 0
          }
        } else if(!is.na(mutationdata[variants[i],]$qflag) & is.finite(p)){
          if(p>=0.05){#if already flagged, but we find a library without a difference between variant and ref scores, remove flag
            mutationdata[variants[i],]$qflag = 0
          }
        }
      } else {
        qualpvalue = 'NA'
      }
      if (length(data[which(data$V3=='ref' & data$V4==samples[j] & data$V5==variants[i]),]$V2)>1 & length(data[which(data$V3=='alt' & data$V4==samples[j] & data$V5==variants[i]),]$V2)>1){
        p = wilcox.test(refmqual, altmqual)$p.value
        mqualpvalue <- paste('mqual p=', round(p,3), sep='')
        if(is.na(mutationdata[variants[i],]$mqflag) & is.finite(p)){
          if(p<0.05){#flag
            mutationdata[variants[i],]$mqflag = 1
          } else {#mark as okay
            mutationdata[variants[i],]$mqflag = 0
          }
        } else if(!is.na(mutationdata[variants[i],]$mqflag) & is.finite(p)){
          if(p>=0.05){#if already flagged, but we find a library without a difference between variant and ref scores, remove flag
            mutationdata[variants[i],]$mqflag = 0
          }
        }
      } else {
        mqualpvalue <- 'NA'
      }
      refCount <- 0
      altCount <- 0
      refQCount <- 0
      altQCount <- 0
      if (length(data_allcounts[which(data_allcounts=='ref' & data_allcounts$V4==samples[j] & data_allcounts$V5==variants[i]),]$freq) > 0){
        refCount <- data_allcounts[which(data_allcounts=='ref' & data_allcounts$V4==samples[j] & data_allcounts$V5==variants[i]),]$freq
      }
      if (length(data_allcounts[which(data_allcounts=='alt' & data_allcounts$V4==samples[j] & data_allcounts$V5==variants[i]),]$freq) > 0){
        altCount <- data_allcounts[which(data_allcounts=='alt' & data_allcounts$V4==samples[j] & data_allcounts$V5==variants[i]),]$freq 
      }
      if (length(data_highqualcounts[which(data_highqualcounts=='ref' & data_highqualcounts$V4==samples[j] & data_highqualcounts$V5==variants[i]),]$freq) > 0){
        refQCount <- data_highqualcounts[which(data_highqualcounts=='ref' & data_highqualcounts$V4==samples[j] & data_highqualcounts$V5==variants[i]),]$freq
      }
      if (length(data_highqualcounts[which(data_highqualcounts=='alt' & data_highqualcounts$V4==samples[j] & data_highqualcounts$V5==variants[i]),]$freq) > 0){
        altQCount <- data_highqualcounts[which(data_highqualcounts=='alt' & data_highqualcounts$V4==samples[j] & data_highqualcounts$V5==variants[i]),]$freq
      }
      allaltcounts = allaltcounts + altQCount
      text(-2,.49,paste(refCount,'/',altCount,' (',refQCount,'/',altQCount,')'), cex=.8, pos=4)
      #compare strands for ref and alt
      a <- 0 #+ ref
      b <- 0 #+ alt
      c <- 0 #- ref
      d <- 0 #- alt
      if (length(data_collapsed_strand[which(data_collapsed_strand$V3=='ref' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='+'),]$freq) > 0){
        a <- data_collapsed_strand[which(data_collapsed_strand$V3=='ref' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='+'),]$freq
      } 
      if (length(data_collapsed_strand[which(data_collapsed_strand$V3=='alt' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='+'),]$freq) > 0){
        b <- data_collapsed_strand[which(data_collapsed_strand$V3=='alt' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='+'),]$freq
      }
      if (length(data_collapsed_strand[which(data_collapsed_strand$V3=='ref' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='-'),]$freq) > 0){
        c <- data_collapsed_strand[which(data_collapsed_strand$V3=='ref' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='-'),]$freq
      } 
      if (length(data_collapsed_strand[which(data_collapsed_strand$V3=='alt' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='-'),]$freq) > 0){
        d <- data_collapsed_strand[which(data_collapsed_strand$V3=='alt' & data_collapsed_strand$V4==samples[j] & data_collapsed_strand$V5==variants[i] & data_collapsed_strand$V6=='-'),]$freq
      }
      if (altQCount>4){
        SBratio <- (b+1)/(d+1)
        if(is.na(mutationdata[variants[i],]$SB)){
          if (SBratio>=2 | SBratio<=.5){
            mutationdata[variants[i],]$SB <- 1
          } else {
            mutationdata[variants[i],]$SB <- 0
          }
        } else if (!is.na(mutationdata[variants[i],]$SB)){
          if (SBratio > .5 & SBratio < 2){
            mutationdata[variants[i],]$SB <- 0
          }
        }
      } else {
        SBratio = NA
      }
      SB <- paste('(+|-) ref:',a,'|',c,'; alt',b,'|',d)
      text(-2,.45,SB, cex=.8, pos=4)
      #compare pair orientation for ref and alt
      a <- 0 #ref F1
      b <- 0 #alt F1
      c <- 0 #ref F2
      d <- 0 #alt F2
      e <- 0 #ref R1
      f <- 0 #alt R1
      g <- 0 #ref R2
      h <- 0 #alt R2
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==1),]$freq) > 0){
        a <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==1),]$freq
      } 
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==1),]$freq) > 0){
        b <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==1),]$freq
      }
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==2),]$freq) > 0){
        c <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==2),]$freq
      } 
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==2),]$freq) > 0){
        d <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='+' & data_collapsed_orientation$V7==2),]$freq
      }
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==1),]$freq) > 0){
        e <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==1),]$freq
      } 
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==1),]$freq) > 0){
        f <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==1),]$freq
      }
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==2),]$freq) > 0){
        g <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='ref' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==2),]$freq
      } 
      if (length(data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==2),]$freq) > 0){
        h <- data_collapsed_orientation[which(data_collapsed_orientation$V3=='alt' & data_collapsed_orientation$V4==samples[j] & data_collapsed_orientation$V5==variants[i] & data_collapsed_orientation$V6=='-' & data_collapsed_orientation$V7==2),]$freq
      }
      refF1R2 <- a+g
      altF1R2 <- b+h
      refF2R1 <- c+e
      altF2R1 <- d+f
      if (altQCount>4){
        OBratio <- (altF2R1+1)/(altF1R2+1)
        if(is.na(mutationdata[variants[i],]$OB)){
          if (OBratio>=2 | OBratio<=.5){
            mutationdata[variants[i],]$OB <- 1
          } else {
            mutationdata[variants[i],]$OB <- 0
          }
        } else if (!is.na(mutationdata[variants[i],]$OB)){
          if (OBratio > .5 & OBratio < 2){
            mutationdata[variants[i],]$OB <- 0
          }
        }
      } else {
        OBratio = NA
      }
      OB <- paste('(1|2) ref:',refF1R2,'|',refF2R1,'; alt',altF1R2,'|',altF2R1)
      text(-2,.41,OB, cex=.8, pos=4)
      text(-2,.37,qualpvalue, cex=.8, pos=4)
      text(-2,.33,mqualpvalue, cex=.8, pos=4)
    }
    if (allaltcounts < 4){
      mutationdata[variants[i],]$fewalt = 1
    }
    dev.off()
  }
}
# OVERALL LIBRARY QUALITY PLOTS

#Overall library quality comparison - quality score
dir.create(paste(patientID,'_qualplots/libraryQuality',sep=''), showWarnings <- FALSE)
pdf(paste(patientID,"_qualplots/libraryQuality/overall_basequality.pdf", sep=""))
ggplot(data, aes(V1, colour = V4)) + geom_density() + theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()
#Overall library quality comparison - mapq score
pdf(paste(patientID,"_qualplots/libraryQuality/overall_mappingquality.pdf", sep=""))+theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(data, aes(V2, colour = V4)) + geom_density()
dev.off()
#Overall library quality for ref vs alt - quality score
pdf(paste(patientID,"_qualplots/libraryQuality/ref_v_alt_basequality.pdf", sep=""))+theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(data, aes(V1, colour = V3)) + geom_density() + facet_wrap('V4')
dev.off()
#Overall library quality for ref vs alt - mapq score
pdf(paste(patientID,"_qualplots/libraryQuality/ref_v_alt_mappingquality.pdf", sep=""))+theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(data, aes(V2, colour = V3)) + geom_density() + facet_wrap('V4')
dev.off()

#VAF PATTERNS and DRIVER PLOTs
dir.create(paste(patientID,'_qualplots/VAFPatterns',sep=''), showWarnings <- FALSE)
# calculate VAFs
vafs <- getVAFs(data, samples, 20, 0)
#Look at estimators of tumor purity and driver mutations (includes drivers of both LG and GBM)
pdf(paste(patientID,"_qualplots/VAFPatterns/commondrivers.pdf", sep=""))
library(dplyr)
drivers <- c('EGFR_','TP53_','PTEN_','NF1_','IDH1_','IDH2_','CIC_','ATRX_','KRAS_','PIK3CA_','PTPRD_','RB1_','MDM2_','PIK3R1_','MSH','MLH','PMS')
banned <- c('IL12RB1_','ZCRB1_','GABRB1_', 'KCNF1_')
cols <- c('#696969',"#33A02C","#E31A1C","#1F78B4","#A6CEE3","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FB9A99","#B2DF8A","#FFFF99","#B15928",'#00F5FF','#556B2F','#7CFC00','#8B0A50','#000080','#FFC125','#D3D3D3')
par(mfrow=c(1,1),mar=c(10,4,2,4))
plot(seq(length(samples)),rep(0,length(samples)), las=2, ylim=c(0,1), col='white', ylab="VAF", xlab='', xaxt="n")
i = 1
j = 1
variantLabels <- c()
notFound <- c()
while (i<=length(drivers)){
  potentialMatchIndexes <- grep(drivers[i],variants)
  bannedMatchIndexes <- c()
  for (entry in banned){
    bannedMatchIndexes <- append(bannedMatchIndexes, grep(entry,variants))
  }
  matches <- variants[!potentialMatchIndexes[potentialMatchIndexes %in% bannedMatchIndexes]]
  if (length(matches)==0){
    notFound <- append(notFound,gsub('_','',drivers[i]))
  }
  for (variantID in matches){
    variantLabels <- append(variantLabels, variantID)
    toPlot <- vafs[variantID,]%>% unlist
    lines(toPlot, las=2, ylim=c(0,.5), col=cols[j], lty=1)
    j = j+1
  }
  i = i+1
}
if (!length(notFound) == length(drivers)){
  axis(1, at=seq(length(samples)),labels=samples, las=2)
  legend(1,1,variantLabels,lty=1,cex=.8,col=cols[1:length(variantLabels)],bty='n')
}
legend('topright',notFound,lty=1,cex=.8,col="white",bty='n', title="Not present:", y.intersp=.8)
dev.off()
#PCA
pdf(paste(patientID,"_qualplots/VAFPatterns/pca.pdf", sep=""))
nrm_count_matrix=as.matrix(vafs)
log=log10(1+nrm_count_matrix)
tlog=t(log)
pca=prcomp(tlog)
summary(pca)
raw <- pca$x[,1:2]
plot(raw[,1], raw[,2],  col='black', xlim=c(min(raw[,1]),max(raw[,1])+.7*max(raw[,1])), pch=20,xlab=paste('PC1 (',toString(round(100*summary(pca)$importance[2,1])),'%)'),ylab=paste('PC2 (',toString(round(100*summary(pca)$importance[2,2])),'%)'),cex=1,cex.lab=1,cex.axis=1)
text(raw[,1], raw[,2], labels=samples, cex=.5, pos=4)
dev.off()
#Frequency across libraries (how many variants are in 1, 2, etc libraries)
pdf(paste(patientID,"_qualplots/VAFPatterns/calldistribution.pdf", sep=""))
detach("package:dplyr", unload=TRUE)
data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff & data_alt$V11=='Y'),]
data_alt_filtered <- data_alt_filtered[c(4,5)]
colnames(data_alt_filtered) <- c('sample','variant')
data_alt_filtered <- unique(data_alt_filtered)
data_alt_filtered <- count(data_alt_filtered, vars = c("variant"))
hist(data_alt_filtered$freq, col='black', freq=FALSE, main='Proportion of variants in N samples',xlab='N',ylab='Proportion of variants')
dev.off()
library(dplyr)

# variant SPECTRA PLOTTING
library(dplyr)
dir.create(paste(patientID,'_qualplots/variantSpectraTest',sep=''), showWarnings <- FALSE)
# Visualization 1: By read
if (length(variants) <= 300){
  detach("package:dplyr", unload=TRUE) #must be done as messes with plyr count function
  data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
  data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff),]
  data_alt_filtered <- data_alt_filtered[c(4,5)]
  colnames(data_alt_filtered) <- c('sample','variant')
  data_alt_filtered$variant <- lapply(data_alt_filtered$variant, as.character)
  for (i in 1:dim(data_alt_filtered)[1]){
    data_alt_filtered[i,2] = extractMuts(as.character(data_alt_filtered[i,2]))
  }
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = c('A>C','T>G','A>G','T>C','A>T','T>A','C>G','G>C','C>A','G>T','C>T','G>A'))
  #plot as counts
  pdf(paste(patientID,"_qualplots/variantSpectra/allReads_counts.pdf", sep=""))
  data_alt_filtered <- count(data_alt_filtered, vars=c("sample", "variant"))
  ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = freq, fill = variant)) + 
    labs(x="Sample",y="Read Counts") + 
    geom_bar(stat="identity") + 
    scale_fill_manual(values = substitutionsColorset) + 
    scale_x_discrete(limits=samples) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
  dev.off()
  #plot as ratios
  pdf(paste(patientID,"_qualplots/variantSpectra/allReads_proportion.pdf", sep=""))
  data_alt_filtered <- ddply(data_alt_filtered, "sample", transform, ratio=getRatio(freq))
  ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = ratio, fill = variant)) + 
    labs(x="Sample",y="Ratio") + 
    geom_bar(stat="identity") + 
    scale_fill_manual(values = substitutionsColorset) + 
    scale_x_discrete(limits=samples) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
  dev.off()
}  
# Visualization 2: As single variantal event
data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
collapse <- 1
data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff),]
data_alt_filtered <- data_alt_filtered[c(4,5)]
colnames(data_alt_filtered) <- c('sample','variant')
data_alt_filtered <- unique(data_alt_filtered)
data_alt_filtered$variant <- lapply(data_alt_filtered$variant, as.character)
for (i in 1:dim(data_alt_filtered)[1]){
  data_alt_filtered[i,2] = extractMuts(as.character(data_alt_filtered[i,2]))
}
if (collapse == 1){
  for (i in 1:dim(data_alt_filtered)[1]){
    data_alt_filtered[i,2] = substitutions_grouping[which(substitutions %in% data_alt_filtered[i,2])]
  }
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = unique(substitutions_grouping))
} else{
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = substitutions) 
}
#plot as counts
pdf(paste(patientID,"_qualplots/variantSpectra/allVariantsPresent_counts.pdf", sep=""))
data_alt_filtered <- count(data_alt_filtered, vars = c("sample", "variant"))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = freq, fill = variant)) + 
  labs(x="Sample",y="Number of variants") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = substitutionsGroupedColorset) + 
  scale_x_discrete(limits=samples) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
  dev.off()
#plot as ratios
pdf(paste(patientID,"_qualplots/variantSpectra/allVariantsPresent_proportion.pdf", sep=""))
data_alt_filtered <- ddply(data_alt_filtered, "sample", transform, ratio=getRatio(freq))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = ratio, fill = variant)) + 
  labs(x="Sample",y="Ratio") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = substitutionsGroupedColorset) + 
  scale_x_discrete(limits=samples) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()

# Visualization 3: By call
data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
collapse <- 1
data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff & data_alt$V11=='Y'),]
data_alt_filtered <- data_alt_filtered[c(4,5)]
colnames(data_alt_filtered) <- c('sample','variant')
data_alt_filtered <- unique(data_alt_filtered)
data_alt_filtered$variant <- lapply(data_alt_filtered$variant, as.character)
for (i in 1:dim(data_alt_filtered)[1]){
  data_alt_filtered[i,2] = extractMuts(as.character(data_alt_filtered[i,2]))
}
if (collapse == 1){
  for (i in 1:dim(data_alt_filtered)[1]){
    data_alt_filtered[i,2] = substitutions_grouping[which(substitutions %in% data_alt_filtered[i,2])]
  }
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = unique(substitutions_grouping))
} else{
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = substitutions) 
}
#plot as counts
pdf(paste(patientID,"_qualplots/variantSpectra/calledVariants_counts.pdf", sep=""))
data_alt_filtered <- count(data_alt_filtered, vars = c("sample", "variant"))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = freq, fill = variant)) + 
  labs(x="Sample",y="Number of variants") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = substitutionsGroupedColorset) + 
  scale_x_discrete(limits=samples) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()
#plot as ratios
pdf(paste(patientID,"_qualplots/variantSpectra/calledVariants_proportion.pdf", sep=""))
data_alt_filtered <- ddply(data_alt_filtered, "sample", transform, ratio=getRatio(freq))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = sample, y = ratio, fill = variant)) + 
  labs(x="Sample",y="Ratio") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = substitutionsGroupedColorset) + 
  scale_x_discrete(limits=samples) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()

# Visualization 4: With trinucleotide content
data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
collapse <- 1
data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff),]
data_alt_filtered = data_alt_filtered[c(4,5,10)]
colnames(data_alt_filtered) <- c('sample','variant','trinucleotide')
data_alt_filtered <- unique(data_alt_filtered)
data_alt_filtered$variant <- lapply(data_alt_filtered$variant, as.character)
for (i in 1:dim(data_alt_filtered)[1]){
  data_alt_filtered[i,2] = extractMuts(as.character(data_alt_filtered[i,2]))
}
if (collapse == 1){
  for (i in 1:dim(data_alt_filtered)[1]){
    data_alt_filtered[i,2] = substitutions_grouping[which(substitutions %in% data_alt_filtered[i,2])]
  }
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = unique(substitutions_grouping))
} else{
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = substitutions) 
}
data_alt_filtered$trinucleotide <-as.factor(as.character(data_alt_filtered$trinucleotide))
data_alt_filtered$trinucleotide <- factor(data_alt_filtered$trinucleotide, levels = trinucleotide_grouping) 
data_alt_filtered$sample <- as.factor(as.character(data_alt_filtered$sample))
data_alt_filtered$sample <- factor(data_alt_filtered$sample, levels = samples)
colourCount = length(unique(data_alt_filtered$trinucleotide))
color = colorRampPalette(brewer.pal(12, "Set3"))(colourCount)
#plot as counts
pdf(paste(patientID,"_qualplots/variantSpectra/allVariantsPresent_trinucleotide_counts.pdf", sep=""))
data_alt_filtered <- count(data_alt_filtered, vars = c("sample", "variant", "trinucleotide"))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = variant, y = freq, fill = trinucleotide)) + 
  labs(x="Substitutions",y="Number of variants") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = trinucleotideColorset) + 
  facet_wrap('sample') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()
#plot as ratios
pdf(paste(patientID,"_qualplots/variantSpectra/allVariantsPresent_trinucleotide_ratios.pdf", sep=""))
data_alt_filtered <- ddply(data_alt_filtered, c("sample","variant"), transform, ratio=getRatio(freq))
ggplot(data = data_alt_filtered[rev(order(data_alt_filtered$variant)),], aes(x = variant, y = ratio, fill = trinucleotide)) + 
  labs(x="Substitutions",y="Proportion with Trinuc") + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = trinucleotideColorset) + 
  facet_wrap('sample') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()

#Visualization 5: Frequency of Substitution Type Across All Samples in Patient
data_alt <- data[which(data$V3=='alt'),]#it doesn't matter if we take mqual or qual, we just want to make sure we only take one
collapse <- 1
data_alt_filtered <- data_alt[which(data_alt$V1>=qual_cutoff & data_alt$V2>=mqual_cutoff & data_alt$V11=='Y'),]
data_alt_filtered <- data_alt_filtered[c(4,5)]
data_alt_filtered <- unique(data_alt_filtered)
data_alt_filtered$variantName <- data_alt_filtered$V5
colnames(data_alt_filtered) <- c('sample','variant','variantName')
data_alt_filtered$variant <- lapply(data_alt_filtered$variant, as.character)
for (i in 1:dim(data_alt_filtered)[1]){
  data_alt_filtered[i,2] = extractMuts(as.character(data_alt_filtered[i,2]))
}
if (collapse == 1){
  for (i in 1:dim(data_alt_filtered)[1]){
    data_alt_filtered[i,2] = substitutions_grouping[which(substitutions %in% data_alt_filtered[i,2])]
  }
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = unique(substitutions_grouping))
} else{
  data_alt_filtered$variant <- as.factor(as.character(data_alt_filtered$variant))
  data_alt_filtered$variant <- factor(data_alt_filtered$variant, levels = substitutions) 
}
data_alt_filtered <- data_alt_filtered[c(2,3)]
#plot as counts
pdf(paste(patientID,"_qualplots/variantSpectra/typefrequency.pdf", sep=""))
data_alt_filtered <- count(data_alt_filtered, vars = c("variant", "variantName"))
ggplot(data_alt_filtered, aes(freq, colour = variant)) + 
  geom_density() +
  scale_fill_manual(values = substitutionsGroupedColorset) +
  theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20), panel.background = element_rect(fill = 'white', colour = 'black'))
dev.off()
