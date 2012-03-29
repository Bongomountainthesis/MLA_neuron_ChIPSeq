
options(stringsAsFactors = FALSE);

library(IRanges)
library(ShortRead)
library(Rsamtools)
library(ChIPpeakAnno)

#region to look up/downstream of tss
topGenes <- 1:50
upstream <- 5000
downstream <- 15000

#region to count overlaps from - should be same as tss regions
minbp <--5000
maxbp <-15000
bin.size <-2 
start <- seq(from=minbp, to=maxbp-bin.size, by=bin.size)
end   <- seq(from=minbp+bin.size-1, to=maxbp, by=bin.size)


####################################################################
##MLA2 K4me3

###need to use ScanBamParam to pull "pos" from bam file
bamFile <- "MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups.bam"

### need to make rd object from bed file

bed <- read.csv(file="all_overlapping_peaks_MLA_H3K4me3.csv")
bed <- bed[,c(6,2,3,4,13,14,15,17,18,11)]

##take top 100 peaks

colnames(bed) <- c("Peak_Name","Peak_Chr","Peak_Start","Peak_End","Tss_Chr","Tss_Start","Tss_End","Tss_Name","Strand","Peak_FDR")

##need to add "chr" to columns...

bed[,"Peak_Chr"] <- paste("chr",bed[,"Peak_Chr"],sep = "")
bed[,"Tss_Chr"] <- paste("chr",bed[,"Tss_Chr"],sep = "")

##need to remove nonunique reads - not entirely sure why they're still in there??

bed.o <- bed[order(bed[,"Peak_FDR"],decreasing=FALSE),]
bed.top <- bed.o[topGenes,]

##take TSS Start and add up/downstream to it

##split into pos/neg strand - for the time being just look at positive strand as it'll cancel profile out otherwise

pos_strand <- bed.top[which(bed.top[,"Strand"] =="1"),]
neg_strand <- bed.top[which(bed.top[,"Strand"] =="-1"),]




############################ TAKE POSITIVE STRANDED WINDOWS FIRST #########

## add/subtract window size
tss_window_start_pos <- pos_strand[,"Tss_Start"]-upstream
tss_window_end_pos <- pos_strand[,"Tss_Start"]+downstream
chr_pos <- pos_strand[,"Peak_Chr"]

# make dataframe of pos tss window
pos_window <- as.data.frame(cbind(chr_pos,tss_window_start_pos,tss_window_end_pos))

pos_tss_window <- cbind(pos_strand[,"Tss_Start"],pos_window)
colnames(pos_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

pos_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(pos_tss_window$Start),
			end = as.numeric(pos_tss_window$End)
			),
			space = pos_tss_window$Chr
		)

what <- c("rname","pos") 
pos_param <- ScanBamParam(which=pos_tss_window.rd, what=what)
pos_bam <- scanBam(bamFile, param=pos_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

pos_relativise.peaks<-function(i){
   cbind(pos_bam[[i]]$pos-pos_tss_window[i,"Tss_Start"],(pos_bam[[i]]$pos+35)-pos_tss_window[i,"Tss_Start"])
}

pos_reads.matrix<-sapply(1:length(pos_bam), pos_relativise.peaks)
pos_reads.matrix<-do.call(rbind,pos_reads.matrix)



############################ TAKE NEGATIVE STRANDED WINDOWS NEXT #########	

tss_window_end_neg <- neg_strand[,"Tss_Start"]+upstream
tss_window_start_neg <- neg_strand[,"Tss_Start"]-downstream
chr_neg <- neg_strand[,"Peak_Chr"]

# make dataframe of neg tss window
neg_window <- as.data.frame(cbind(chr_neg,tss_window_start_neg,tss_window_end_neg))

neg_tss_window <- cbind(neg_strand[,"Tss_Start"],neg_window)
colnames(neg_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

neg_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(neg_tss_window$Start),
			end = as.numeric(neg_tss_window$End)
			),
			space = neg_tss_window$Chr
		)

what <- c("rname","pos") 
neg_param <- ScanBamParam(which=neg_tss_window.rd, what=what)
neg_bam <- scanBam(bamFile, param=neg_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

neg_relativise.peaks<-function(i){
   cbind(neg_bam[[i]]$pos-neg_tss_window[i,"Tss_Start"],(neg_bam[[i]]$pos+35)-neg_tss_window[i,"Tss_Start"])
}

###for some reason when only one peak exists, sapply builds a matrix instead of a list, so fudge this...
neg_reads.matrix<-sapply(1:length(neg_bam), neg_relativise.peaks)
neg_reads.matrix<-do.call(rbind,neg_reads.matrix)

##switch strands - multiply by -1 and change orientation

neg_reads.matrix <- neg_reads.matrix * -1

neg_reads.matrix <- neg_reads.matrix[,c(2,1)]

################################ cat positive and negative strands back together

##need to check that there were reads in each before combining, as something wierd happens when there are no reads in ranged data....
#this appears to be a massive fudge too - see if matrices are NA, if not combine...

if(!is.na(pos_reads.matrix[1,1]) && !is.na(neg_reads.matrix[1,1])) reads.matrix <- as.data.frame(rbind(pos_reads.matrix, neg_reads.matrix)) else reads.matrix <- as.data.frame(pos_reads.matrix)

colnames(reads.matrix) <- c("Start","End")

###make reads RD and overlap to v_tss

reads <- RangedData(ranges = IRanges(
			start = reads.matrix[,"Start"],
			end = reads.matrix[,"End"]
			)
		)

###need to take peaks only over tss - ie go back to peak list and pull out those bam files...

#count overlaps of virtual tss 50bp regions with peaks

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )


overlap <- countOverlaps(v_tss,
                              reads,
                              type = "within"
                              )[[1]]

overlap_average <- overlap/length(c(pos_bam,neg_bam))

#overlap_average <- overlap/nrow(reads.matrix)
overlap_to_tss_MLA_k4 <- cbind(start, overlap_average)

#######################################################################################################
##MLA2 K27me3

###need to use ScanBamParam to pull "pos" from bam file
bamFile <- "MLA_NS_H3K27me3_CMN056_s_4_export_sorted_nodups.bam"

### need to make rd object from bed file

bed <- read.csv(file="all_overlapping_peaks_MLA_H3K27me3.csv")
bed <- bed[,c(6,2,3,4,13,14,15,17,18,11)]

##take top 100 peaks

colnames(bed) <- c("Peak_Name","Peak_Chr","Peak_Start","Peak_End","Tss_Chr","Tss_Start","Tss_End","Tss_Name","Strand","Peak_FDR")

##need to add "chr" to columns...

bed[,"Peak_Chr"] <- paste("chr",bed[,"Peak_Chr"],sep = "")
bed[,"Tss_Chr"] <- paste("chr",bed[,"Tss_Chr"],sep = "")

##need to remove nonunique reads - not entirely sure why they're still in there??

bed.o <- bed[order(bed[,"Peak_FDR"],decreasing=FALSE),]
bed.top <- bed.o[topGenes,]

##take TSS Start and add up/downstream to it

##split into pos/neg strand - for the time being just look at positive strand as it'll cancel profile out otherwise

pos_strand <- bed.top[which(bed.top[,"Strand"] =="1"),]
neg_strand <- bed.top[which(bed.top[,"Strand"] =="-1"),]




############################ TAKE POSITIVE STRANDED WINDOWS FIRST #########

## add/subtract window size
tss_window_start_pos <- pos_strand[,"Tss_Start"]-upstream
tss_window_end_pos <- pos_strand[,"Tss_Start"]+downstream
chr_pos <- pos_strand[,"Peak_Chr"]

# make dataframe of pos tss window
pos_window <- as.data.frame(cbind(chr_pos,tss_window_start_pos,tss_window_end_pos))

pos_tss_window <- cbind(pos_strand[,"Tss_Start"],pos_window)
colnames(pos_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

pos_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(pos_tss_window$Start),
			end = as.numeric(pos_tss_window$End)
			),
			space = pos_tss_window$Chr
		)

what <- c("rname","pos") 
pos_param <- ScanBamParam(which=pos_tss_window.rd, what=what)
pos_bam <- scanBam(bamFile, param=pos_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

pos_relativise.peaks<-function(i){
   cbind(pos_bam[[i]]$pos-pos_tss_window[i,"Tss_Start"],(pos_bam[[i]]$pos+35)-pos_tss_window[i,"Tss_Start"])
}

pos_reads.matrix<-sapply(1:length(pos_bam), pos_relativise.peaks)
pos_reads.matrix<-do.call(rbind,pos_reads.matrix)



############################ TAKE NEGATIVE STRANDED WINDOWS NEXT #########	

tss_window_end_neg <- neg_strand[,"Tss_Start"]+upstream
tss_window_start_neg <- neg_strand[,"Tss_Start"]-downstream
chr_neg <- neg_strand[,"Peak_Chr"]

# make dataframe of neg tss window
neg_window <- as.data.frame(cbind(chr_neg,tss_window_start_neg,tss_window_end_neg))

neg_tss_window <- cbind(neg_strand[,"Tss_Start"],neg_window)
colnames(neg_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

neg_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(neg_tss_window$Start),
			end = as.numeric(neg_tss_window$End)
			),
			space = neg_tss_window$Chr
		)

what <- c("rname","pos") 
neg_param <- ScanBamParam(which=neg_tss_window.rd, what=what)
neg_bam <- scanBam(bamFile, param=neg_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

neg_relativise.peaks<-function(i){
   cbind(neg_bam[[i]]$pos-neg_tss_window[i,"Tss_Start"],(neg_bam[[i]]$pos+35)-neg_tss_window[i,"Tss_Start"])
}

neg_reads.matrix<-sapply(1:length(neg_bam), neg_relativise.peaks)
neg_reads.matrix<-do.call(rbind,neg_reads.matrix)

##switch strands - multiply by -1 and change orientation

neg_reads.matrix <- neg_reads.matrix * -1

neg_reads.matrix <- neg_reads.matrix[,c(2,1)]

################################ cat positive and negative strands back together


##need to check that there were reads in each before combining, as something wierd happens when there are no reads in ranged data....
#this appears to be a massive fudge too - see if matrices are NA, if not combine...

if(!is.na(pos_reads.matrix[1,1]) && !is.na(neg_reads.matrix[1,1])) reads.matrix <- as.data.frame(rbind(pos_reads.matrix, neg_reads.matrix)) else reads.matrix <- as.data.frame(pos_reads.matrix)

colnames(reads.matrix) <- c("Start","End")

###make reads RD and overlap to v_tss

reads <- RangedData(ranges = IRanges(
			start = reads.matrix[,"Start"],
			end = reads.matrix[,"End"]
			)
		)

###need to take peaks only over tss - ie go back to peak list and pull out those bam files...

#count overlaps of virtual tss 50bp regions with peaks

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )


overlap <- countOverlaps(v_tss,
                              reads,
                              type = "within"
                              )[[1]]

overlap_average <- overlap/length(c(pos_bam,neg_bam))

#overlap_average <- overlap/nrow(reads.matrix)
overlap_to_tss_MLA_k27 <- cbind(start, overlap_average)

#########################################################################################################
##Neuron K4me3


###need to use ScanBamParam to pull "pos" from bam file
bamFile <- "Neuron_H3K4me3_CMN058_s_6_export_sorted_nodups.bam"

### need to make rd object from bed file

bed <- read.csv(file="all_overlapping_peaks_Neuron_H3K4me3.csv")
bed <- bed[,c(6,2,3,4,13,14,15,17,18,11)]

##take top 100 peaks

colnames(bed) <- c("Peak_Name","Peak_Chr","Peak_Start","Peak_End","Tss_Chr","Tss_Start","Tss_End","Tss_Name","Strand","Peak_FDR")

##need to add "chr" to columns...

bed[,"Peak_Chr"] <- paste("chr",bed[,"Peak_Chr"],sep = "")
bed[,"Tss_Chr"] <- paste("chr",bed[,"Tss_Chr"],sep = "")

##need to remove nonunique reads - not entirely sure why they're still in there??

bed.o <- bed[order(bed[,"Peak_FDR"],decreasing=FALSE),]

#for some reason got a load of crap peaks at the top (FDR of 0??)- remove these
#bed.o0 <- bed[which(bed.o[,"Peak_FDR"] > 0),]
bed.top <- bed.o[topGenes,]

##take TSS Start and add up/downstream to it

##split into pos/neg strand - for the time being just look at positive strand as it'll cancel profile out otherwise

pos_strand <- bed.top[which(bed.top[,"Strand"] =="1"),]
neg_strand <- bed.top[which(bed.top[,"Strand"] =="-1"),]




############################ TAKE POSITIVE STRANDED WINDOWS FIRST #########

## add/subtract window size
tss_window_start_pos <- pos_strand[,"Tss_Start"]-upstream
tss_window_end_pos <- pos_strand[,"Tss_Start"]+downstream
chr_pos <- pos_strand[,"Peak_Chr"]

# make dataframe of pos tss window
pos_window <- as.data.frame(cbind(chr_pos,tss_window_start_pos,tss_window_end_pos))

pos_tss_window <- cbind(pos_strand[,"Tss_Start"],pos_window)
colnames(pos_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

pos_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(pos_tss_window$Start),
			end = as.numeric(pos_tss_window$End)
			),
			space = pos_tss_window$Chr
		)

what <- c("rname","pos") 
pos_param <- ScanBamParam(which=pos_tss_window.rd, what=what)
pos_bam <- scanBam(bamFile, param=pos_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

pos_relativise.peaks<-function(i){
   cbind(pos_bam[[i]]$pos-pos_tss_window[i,"Tss_Start"],(pos_bam[[i]]$pos+35)-pos_tss_window[i,"Tss_Start"])
}

pos_reads.matrix<-sapply(1:length(pos_bam), pos_relativise.peaks)
pos_reads.matrix<-do.call(rbind,pos_reads.matrix)



############################ TAKE NEGATIVE STRANDED WINDOWS NEXT #########	

tss_window_end_neg <- neg_strand[,"Tss_Start"]+upstream
tss_window_start_neg <- neg_strand[,"Tss_Start"]-downstream
chr_neg <- neg_strand[,"Peak_Chr"]

# make dataframe of neg tss window
neg_window <- as.data.frame(cbind(chr_neg,tss_window_start_neg,tss_window_end_neg))

neg_tss_window <- cbind(neg_strand[,"Tss_Start"],neg_window)
colnames(neg_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

neg_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(neg_tss_window$Start),
			end = as.numeric(neg_tss_window$End)
			),
			space = neg_tss_window$Chr
		)

what <- c("rname","pos") 
neg_param <- ScanBamParam(which=neg_tss_window.rd, what=what)
neg_bam <- scanBam(bamFile, param=neg_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

neg_relativise.peaks<-function(i){
   cbind(neg_bam[[i]]$pos-neg_tss_window[i,"Tss_Start"],(neg_bam[[i]]$pos+35)-neg_tss_window[i,"Tss_Start"])
}

###for some reason when only one peak exists, sapply builds a matrix instead of a list, so fudge this...
neg_reads.matrix<-sapply(1:length(neg_bam), neg_relativise.peaks)
neg_reads.matrix<-do.call(rbind,neg_reads.matrix)

##switch strands - multiply by -1 and change orientation

neg_reads.matrix <- neg_reads.matrix * -1

neg_reads.matrix <- neg_reads.matrix[,c(2,1)]

################################ cat positive and negative strands back together

##need to check that there were reads in each before combining, as something wierd happens when there are no reads in ranged data....
#this appears to be a massive fudge too - see if matrices are NA, if not combine...

if(!is.na(pos_reads.matrix[1,1]) && !is.na(neg_reads.matrix[1,1])) reads.matrix <- as.data.frame(rbind(pos_reads.matrix, neg_reads.matrix)) else reads.matrix <- as.data.frame(pos_reads.matrix)

colnames(reads.matrix) <- c("Start","End")

###make reads RD and overlap to v_tss

reads <- RangedData(ranges = IRanges(
			start = reads.matrix[,"Start"],
			end = reads.matrix[,"End"]
			)
		)

###need to take peaks only over tss - ie go back to peak list and pull out those bam files...

#count overlaps of virtual tss 50bp regions with peaks

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )


overlap <- countOverlaps(v_tss,
                              reads,
                              type = "within"
                              )[[1]]

#overlap_average <- overlap/length(c(pos_bam,neg_bam))

#overlap_average <- overlap/nrow(reads.matrix)
overlap_to_tss_neuron_k4 <- cbind(start, overlap)


###############################################################################################################
##Neuron K27me3

###need to use ScanBamParam to pull "pos" from bam file
bamFile <- "Neuron2_H3K27me3_CMN059_s_7_export_sorted_nodups.bam"

### need to make rd object from bed file

bed <- read.csv(file="all_overlapping_peaks_Neuron_H3K27me3.csv")
bed <- bed[,c(6,2,3,4,13,14,15,17,18,11)]

##take top 100 peaks

colnames(bed) <- c("Peak_Name","Peak_Chr","Peak_Start","Peak_End","Tss_Chr","Tss_Start","Tss_End","Tss_Name","Strand","Peak_FDR")

##need to add "chr" to columns...

bed[,"Peak_Chr"] <- paste("chr",bed[,"Peak_Chr"],sep = "")
bed[,"Tss_Chr"] <- paste("chr",bed[,"Tss_Chr"],sep = "")

##need to remove nonunique reads - not entirely sure why they're still in there??

bed.o <- bed[order(bed[,"Peak_FDR"],decreasing=FALSE),]
bed.top <- bed.o[topGenes,]

##take TSS Start and add up/downstream to it

##split into pos/neg strand - for the time being just look at positive strand as it'll cancel profile out otherwise

pos_strand <- bed.top[which(bed.top[,"Strand"] =="1"),]
neg_strand <- bed.top[which(bed.top[,"Strand"] =="-1"),]




############################ TAKE POSITIVE STRANDED WINDOWS FIRST #########

## add/subtract window size
tss_window_start_pos <- pos_strand[,"Tss_Start"]-upstream
tss_window_end_pos <- pos_strand[,"Tss_Start"]+downstream
chr_pos <- pos_strand[,"Peak_Chr"]

# make dataframe of pos tss window
pos_window <- as.data.frame(cbind(chr_pos,tss_window_start_pos,tss_window_end_pos))

pos_tss_window <- cbind(pos_strand[,"Tss_Start"],pos_window)
colnames(pos_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

pos_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(pos_tss_window$Start),
			end = as.numeric(pos_tss_window$End)
			),
			space = pos_tss_window$Chr
		)

what <- c("rname","pos") 
pos_param <- ScanBamParam(which=pos_tss_window.rd, what=what)
pos_bam <- scanBam(bamFile, param=pos_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

pos_relativise.peaks<-function(i){
   cbind(pos_bam[[i]]$pos-pos_tss_window[i,"Tss_Start"],(pos_bam[[i]]$pos+35)-pos_tss_window[i,"Tss_Start"])
}

pos_reads.matrix<-sapply(1:length(pos_bam), pos_relativise.peaks)
pos_reads.matrix<-do.call(rbind,pos_reads.matrix)



############################ TAKE NEGATIVE STRANDED WINDOWS NEXT #########	

tss_window_end_neg <- neg_strand[,"Tss_Start"]+upstream
tss_window_start_neg <- neg_strand[,"Tss_Start"]-downstream
chr_neg <- neg_strand[,"Peak_Chr"]

# make dataframe of neg tss window
neg_window <- as.data.frame(cbind(chr_neg,tss_window_start_neg,tss_window_end_neg))

neg_tss_window <- cbind(neg_strand[,"Tss_Start"],neg_window)
colnames(neg_tss_window) <- c("Tss_Start","Chr","Start","End")

##then make ranged data object from the tss_window

neg_tss_window.rd <- RangedData(ranges = IRanges(
			start = as.numeric(neg_tss_window$Start),
			end = as.numeric(neg_tss_window$End)
			),
			space = neg_tss_window$Chr
		)

what <- c("rname","pos") 
neg_param <- ScanBamParam(which=neg_tss_window.rd, what=what)
neg_bam <- scanBam(bamFile, param=neg_param)

##dont need this, but just in case want read positions...
#bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

##make read start position relative to TSS start position

neg_relativise.peaks<-function(i){
   cbind(neg_bam[[i]]$pos-neg_tss_window[i,"Tss_Start"],(neg_bam[[i]]$pos+35)-neg_tss_window[i,"Tss_Start"])
}

neg_reads.matrix<-sapply(1:length(neg_bam), neg_relativise.peaks)
neg_reads.matrix<-do.call(rbind,neg_reads.matrix)

##switch strands - multiply by -1 and change orientation

neg_reads.matrix <- neg_reads.matrix * -1

neg_reads.matrix <- neg_reads.matrix[,c(2,1)]

################################ cat positive and negative strands back together


##need to check that there were reads in each before combining, as something wierd happens when there are no reads in ranged data....
#this appears to be a massive fudge too - see if matrices are NA, if not combine...

if(!is.na(pos_reads.matrix[1,1]) && !is.na(neg_reads.matrix[1,1])) reads.matrix <- as.data.frame(rbind(pos_reads.matrix, neg_reads.matrix)) else reads.matrix <- as.data.frame(pos_reads.matrix)

colnames(reads.matrix) <- c("Start","End")

###make reads RD and overlap to v_tss

reads <- RangedData(ranges = IRanges(
			start = reads.matrix[,"Start"],
			end = reads.matrix[,"End"]
			)
		)

###need to take peaks only over tss - ie go back to peak list and pull out those bam files...

#count overlaps of virtual tss 50bp regions with peaks

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )


overlap <- countOverlaps(v_tss,
                              reads,
                              type = "within"
                              )[[1]]

#overlap_average <- overlap/length(c(pos_bam,neg_bam))

#overlap_average <- overlap/nrow(reads.matrix)
overlap_to_tss_neuron_k27 <- cbind(start, overlap)

##########################################################################################################

legend(75000, 120000, c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)


##draw phat plot of all peak shapes


postscript("peak_shape_all.ps", paper="special",width=10,height=10,horizontal = FALSE)
plot(overlap_to_tss_MLA_k4[,"start"], 
     overlap_to_tss_MLA_k4[,"overlap_average"], 
     col="navy", 
     main="Distribution of histone enrichment peaks around the TSS", 
     xlab= "Distance to TSS (bp)", 
     ylab="Relative number of reads",
     type="l",
     lwd = 3)

lines(overlap_to_tss_MLA_k27[,"start"],
      overlap_to_tss_MLA_k27[,"overlap_average"],
      col="firebrick1",
      lwd = 3)

legend(9000,10.5, c("H3K4me3","H3K27me3"), col = c("firebrick1","navy"), lty=1, lwd = 5)
#lines(overlap_to_tss_neuron_k4[,"start"],
#      overlap_to_tss_neuron_k4[,"overlap"],
#      col="blue",
#      type="l"
#      )
#lines(overlap_to_tss_neuron_k27[,"start"],
#      overlap_to_tss_neuron_k27[,"overlap"],
#      col="brown",
#      type="l"
#      )

dev.off()



plot(overlap_to_tss_neuron_k4[,"start"], 
     overlap_to_tss_neuron_k4[,"overlap"], 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="Number of Reads",
     type="l"
)
