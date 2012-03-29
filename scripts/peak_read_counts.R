#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
###
#
# For each region in the BED files, fetches the read count in that regions from the BAM file.
#
###

##bed files

biv_ns <- read.csv(file = "/mnt/data/bivalent_NS_genes.csv")
biv_n <- read.csv(file = "/mnt/data/bivalent_N_genes.csv")

##need to remove duplicate probes names otherwise it'll all be abit screwed...

biv_ns <- biv_ns[order(abs(biv_ns[,"logFC"]),decreasing = TRUE),]
biv_ns <- biv_ns[!duplicated(biv_ns[,"NS_K27"]),]

biv_n <- biv_n[order(abs(biv_n[,"logFC"]),decreasing = TRUE),]
biv_n <- biv_n[!duplicated(biv_n[,"N_K27"]),]



####### just deal with NSC bivalent genes for now - take NS_K27 column and split it out

biv_ns_peaks <- biv_ns[,"NS_K27"]

biv_ns_chr <- gsub(":.*$","",biv_ns_peaks)

biv_ns_start <- gsub("^.*:","",biv_ns_peaks)
biv_ns_starter <- gsub("-.*$","",biv_ns_start)

biv_ns_end <- gsub("^.*-","",biv_ns_peaks)

##cbind back together

biv_ns_peaks_sep <- as.data.frame(cbind(biv_ns_chr, biv_ns_starter, biv_ns_end))

colnames(biv_ns_peaks_sep) <- c("Peak_chr","Peak_start","Peak_end")

rownames(biv_ns_peaks_sep) <- biv_ns_peaks

##remove 1 from start as bed starts at 0

biv_ns_peaks_sep[,"Peak_start"] <- as.numeric(biv_ns_peaks_sep[,"Peak_start"]) - 1

write.table(biv_ns_peaks_sep, file = "/mnt/data/NS_bivalent_peaks_for_deseq.bed",row.names = F, col.names = F, quote = F)

#######maybe take N_K27 too - bit more of an arse. Take Neuron K27 peak if there, if not fill in with NSC one
#biv_n_res <- biv_ns
#
#for(i in 1:length(biv_n_res[,"NS_K4"])){
# 	if (is.na(biv_n_res[i,"N_K27"])){
#		biv_n_res[i,"N_K27"] <- biv_ns[i,"NS_K27"]
#	} else {
#		biv_n_res[i,"N_K27"] <- biv_n_res[i,"N_K27"]
#	}
#}

## split up as before

biv_n_peaks <- biv_n[,"N_K27"]

biv_n_chr <- gsub(":.*$","",biv_n_peaks)

biv_n_start <- gsub("^.*:","",biv_n_peaks)
biv_n_starter <- gsub("-.*$","",biv_n_start)

biv_n_end <- gsub("^.*-","",biv_n_peaks)

##cbind back together

biv_n_peaks_sep <- as.data.frame(cbind(biv_n_chr, biv_n_starter, biv_n_end))

colnames(biv_n_peaks_sep) <- c("Peak_chr","Peak_start","Peak_end")

rownames(biv_n_peaks_sep) <- biv_n_peaks

##remove 1 from start as bed starts at 0

biv_n_peaks_sep[,"Peak_start"] <- as.numeric(biv_n_peaks_sep[,"Peak_start"]) - 1

write.table(biv_n_peaks_sep, file = "/mnt/data/N_bivalent_peaks_for_deseq.bed",row.names = F, col.names = F, quote = F)

############just in case, should calculate middle region of peak whilst here - blatantly will need that bit instead

biv_ns_peaks_mid <-biv_ns_peaks_sep

window <-3000 

for(i in 1:length(biv_ns_peaks_mid[,1])){
	biv_ns_peaks_mid[i,4] <- (as.numeric(biv_ns_peaks_mid[i,"Peak_end"]) - as.numeric(biv_ns_peaks_mid[i,"Peak_start"])) / 2
	biv_ns_peaks_mid[i,"Peak_start"] <- (as.numeric(biv_ns_peaks_mid[i,"Peak_start"]) + as.numeric(biv_ns_peaks_mid[i,4])) - 250
	biv_ns_peaks_mid[i,"Peak_end"] <- (as.numeric(biv_ns_peaks_mid[i,"Peak_end"]) - as.numeric(biv_ns_peaks_mid[i,4])) + 250
}

biv_n_peaks_mid <-biv_n_peaks_sep

for(i in 1:length(biv_n_peaks_mid[,1])){
	biv_n_peaks_mid[i,4] <- (as.numeric(biv_n_peaks_mid[i,"Peak_end"]) - as.numeric(biv_n_peaks_mid[i,"Peak_start"])) / 2
	biv_n_peaks_mid[i,"Peak_start"] <- (as.numeric(biv_n_peaks_mid[i,"Peak_start"]) + as.numeric(biv_n_peaks_mid[i,4])) - 250
	biv_n_peaks_mid[i,"Peak_end"] <- (as.numeric(biv_n_peaks_mid[i,"Peak_end"]) - as.numeric(biv_n_peaks_mid[i,4])) + 250
}

write.table(biv_ns_peaks_mid, file = "/mnt/data/NS_bivalent_peaks_middle_for_deseq.bed",row.names = F, col.names = F, quote = F)
write.table(biv_n_peaks_mid, file = "/mnt/data/N_bivalent_peaks_middle_for_deseq.bed",row.names = F, col.names = F, quote = F)


# you will need:
# The BAM files for each of the samples
# A BED file containing the peak regions for each of the samples


#!/usr/local/bin/Rscript peak_read_counts.R /path/to/outdir N_threads /path/to/bam1 /path/to/bed1 /path/to/bam2 /path/to/bed2 ...

options(stringsAsFactors = FALSE);


args <- commandArgs(trailingOnly=TRUE)

#for testing
#args<-c("/mnt/data/",4, "/mnt/astro/BAM/IP.bam" ,"/mnt/astro/Macs/NA_peaks.bed", "/mnt/esc/chip_export_sorted_nodups.bam","/mnt/esc/macs_300_1.0e-05/EscChIPseqREST_peaks.bed")
#args <- c("/mnt/data/",4,"/mnt/TC/CMN066_s_8_export_sorted_nodups.bam", "/mnt/TC/macs_300_1.0e-05/Mash1_TC_peaks.bed", "/mnt/SC/chip_export_sorted_nodups.bam", "/mnt/SC/macs_300_1.0e-05/Mash1_SC_peaks.bed")

args<-c("/mnt/data",4, "MLA_NS_H3K27me3_CMN056_s_4_export_sorted_nodups.bam","/mnt/data/NS_bivalent_peaks_for_deseq.bed", "/mnt/data/Neuron2_H3K27me3_CMN059_s_7_export_sorted_nodups.bam","/mnt/data/N_bivalent_peaks_for_deseq.bed")
#args<-c("/mnt/data",4, "/mnt/data/Neuron2_H3K27me3_CMN059_s_7_export_sorted_nodups.bam","/mnt/data/N_bivalent_peaks_middle_for_deseq.bed", "MLA_NS_H3K27me3_CMN056_s_4_export_sorted_nodups.bam","/mnt/data/NS_bivalent_peaks_middle_for_deseq.bed")

resdir <- args[1]
threads <- args[2]
args <- args[-1:-2]

inds <- 1:length(args)
bam.files <- args[which(inds%%2!=0)]
bed.files <- args[which(inds%%2==0)]


library(IRanges)
library(ShortRead)
library(snow)
library(baySeq)
library(DESeq)
library(rtracklayer)
library(Rsamtools)

# For ChIPseq data we aren't dealing with that many locations, so probably
# we're good using just the number of cores on the AWS machine. 
#if(is.null(threads) || threads==1){
#  cl <- NULL
#}else{
#  cl <- makeCluster(threads,"SOCK")
#}


# Read in the BED files as RangedData.
# Retrieve the read data for those regions
# Build the counts table
# This seems to take about an hour for 4K regions.

beds <- list()
counts <- NULL
seglens <- NULL
for(i in 1:length(bed.files)){
  beds[[i]] <- import(bed.files[i])

  #get the count data for these ranges from each bam file
  bam.counts <- list()
  for(j in 1:length(bam.files)){
    what <- c("qname") 
    param <- ScanBamParam(which=beds[[i]], what=what)
    bam <- scanBam(bam.files[[j]], param=param)
    bam.counts[[j]] <- sapply(bam, function(x){length(x$qname)})
  }

  beds[[i]] <- as.data.frame(beds[[i]])
  nms <- paste(beds[[i]][,"space"], paste(beds[[i]][,"start"], beds[[i]][,"end"], sep="-"), sep=":")
  these <- do.call(cbind, bam.counts)
  col.nms <- sub("-","", gsub("/","-", bam.files))
  colnames(these) <- col.nms
  these <- data.frame(sample=bed.files[i], these)
  if(is.null(counts)){
    counts <- these
  }else{
    counts <- rbind(counts, these)
  }
  if(is.null(seglens)){
    seglens <- beds[[i]][,"end"]-beds[[i]][,"start"]+1
  }else{
    seglens <- c(seglens,beds[[i]][,"end"]-beds[[i]][,"start"]+1)
  }
}

save(counts, file=paste(resdir,"counts.RData", sep="/"))



# We'll need the library sizes later, so get them from the summary files
# Note that these should be generated using the samtools flagstat program
summaries <- gsub('.bam', '.summary', bam.files )
libsizes <- sapply(summaries, function(x){
  l <- readLines(x,1)
  l <- sub("\\s+.*","",l, perl=T)
})


save(libsizes, file=paste(resdir,"libsizes.RData", sep="/"))
