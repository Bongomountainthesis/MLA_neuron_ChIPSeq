#!/usr/local/bin/Rscript


###load limma results file and nearest or overlapping feature data and merge them together

limma.filename <- "/space/matt/mla_neurons_lumixpn/results/limma_results.csv"
NS.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K4me3/all_overlapping_peaks.csv"
NS.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K27me3/all_overlapping_peaks.csv"
N.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K4me3/all_overlapping_peaks.csv"
N.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K27me3/all_overlapping_peaks.csv"

limma <- read.csv(limma.filename)
NS.K4 <- read.csv(NS.K4.filename)
NS.K27 <- read.csv(NS.K27.filename)
N.K4 <- read.csv(N.K4.filename)
N.K27 <- read.csv(N.K27.filename)

##remove duplicate probe-gene matches - keep highest FC difference
limma.o <- limma[order(abs(limma[,"logFC"]), decreasing=TRUE),]
limma.do <- limma.o[!duplicated(limma.o[,"EnsemblID"]),]



##############plot density plot of sizes of peaks

postscript("peak_width.ps", horizontal=FALSE)
plot(density(N.K4[,"rd.width"]), 
             col = "green",
             main = "H3K4me3 and H3K27me3 Modification Peak Widths", 
             xlab = "Peak Width (bp)", 
             ylab = "Number of Peaks",
             xlim = c(0, 80000)
             
          )

lines(col = "red", density(NS.K27[,"rd.width"]))
lines(col = "blue", density(NS.K4[,"rd.width"]))
lines(col = "magenta", density(N.K27[,"rd.width"]))

legend(60000, 5e-04 , c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)

dev.off()


##################plot numbers of peaks over a virtual tss
#calculate distances from peak start/end to tss
###strandedness - work out which rows are positive and negative


library(IRanges)
NS.K4.pos <- which(NS.K4[,"tss.strand"]==1)
NS.K4.neg <- which(NS.K4[,"tss.strand"]==-1)

NS.K4.peakstart <- c((NS.K4[NS.K4.pos,"rd.start"] - NS.K4[NS.K4.pos,"tss.start"]),(NS.K4[NS.K4.neg,"tss.end"] - NS.K4[NS.K4.neg,"rd.end"]))
NS.K4.peakend <- c((NS.K4[NS.K4.pos,"rd.end"] - NS.K4[NS.K4.pos,"tss.start"]),(NS.K4[NS.K4.neg,"tss.end"] - NS.K4[NS.K4.neg,"rd.start"]))
NS.K4.ranges <- cbind(NS.K4.peakstart,NS.K4.peakend)

NS.K4.ir <- RangedData(ranges = IRanges(
                       start = NS.K4.ranges[,1],
                       end = NS.K4.ranges[,2],
                       )
                    )          

NS.K27.pos <- which(NS.K27[,"tss.strand"]==1)
NS.K27.neg <- which(NS.K27[,"tss.strand"]==-1)

NS.K27.peakstart <- c((NS.K27[NS.K27.pos,"rd.start"] - NS.K27[NS.K27.pos,"tss.start"]),(NS.K27[NS.K27.neg,"tss.end"] - NS.K27[NS.K27.neg,"rd.end"]))
NS.K27.peakend <- c((NS.K27[NS.K27.pos,"rd.end"] - NS.K27[NS.K27.pos,"tss.start"]),(NS.K27[NS.K27.neg,"tss.end"] - NS.K27[NS.K27.neg,"rd.start"]))
NS.K27.ranges <- cbind(NS.K27.peakstart,NS.K27.peakend)
NS.K27.ir <- RangedData(ranges = IRanges(
                       start = NS.K27.ranges[,1],
                       end = NS.K27.ranges[,2],
                       )
                    )

N.K4.pos <- which(N.K4[,"tss.strand"]==1)
N.K4.neg <- which(N.K4[,"tss.strand"]==-1)

N.K4.peakstart <- c((N.K4[N.K4.pos,"rd.start"] - N.K4[N.K4.pos,"tss.start"]),(N.K4[N.K4.neg,"tss.end"] - N.K4[N.K4.neg,"rd.end"]))
N.K4.peakend <- c((N.K4[N.K4.pos,"rd.end"] - N.K4[N.K4.pos,"tss.start"]),(N.K4[N.K4.neg,"tss.end"] - N.K4[N.K4.neg,"rd.start"]))
N.K4.ranges <- cbind(N.K4.peakstart,N.K4.peakend)
N.K4.ir <- RangedData(ranges = IRanges(
                       start = N.K4.ranges[,1],
                       end = N.K4.ranges[,2],
                       )
                    )

N.K27.pos <- which(N.K27[,"tss.strand"]==1)
N.K27.neg <- which(N.K27[,"tss.strand"]==-1)

N.K27.peakstart <- c((N.K27[N.K27.pos,"rd.start"] - N.K27[N.K27.pos,"tss.start"]),(N.K27[N.K27.neg,"tss.end"] - N.K27[N.K27.neg,"rd.end"]))
N.K27.peakend <- c((N.K27[N.K27.pos,"rd.end"] - N.K27[N.K27.pos,"tss.start"]),(N.K27[N.K27.neg,"tss.end"] - N.K27[N.K27.neg,"rd.start"]))
N.K27.ranges <- cbind(N.K27.peakstart,N.K27.peakend)
N.K27.ir <- RangedData(ranges = IRanges(
                       start = N.K27.ranges[,1],
                       end = N.K27.ranges[,2],
                       )
                    )

#count overlaps of virtual tss 50bp regions with peaks

minbp <--100000
maxbp <-100000
bin.size <-50 
start <- seq(from=minbp, to=maxbp-bin.size, by=bin.size)
end   <- seq(from=minbp+bin.size-1, to=maxbp, by=bin.size)

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )

NS.K4.overlap <- findOverlaps(v_tss,
                              NS.K4.ir,
                              type = "within"
                              
                             )
NS.K4.overlap <- as.matrix(NS.K4.overlap[[1]])

NS.K27.overlap <- findOverlaps(v_tss, 
                              NS.K27.ir, 
                              type = "within", 
                             )
NS.K27.overlap <- as.matrix(NS.K27.overlap[[1]])

N.K4.overlap <- findOverlaps(v_tss, 
                              N.K4.ir, 
                              type = "within", 
                             )
N.K4.overlap <- as.matrix(N.K4.overlap[[1]])

N.K27.overlap <- findOverlaps(v_tss, 
                              N.K27.ir, 
                              type = "within", 
                             )
N.K27.overlap <- as.matrix(N.K27.overlap[[1]])


##replace position in original (e.g. NS.K4) dataframe with "rd.values.FoldEnrichment"
NS.K4.fe <- cbind(NS.K4.overlap, NS.K4[NS.K4.overlap[,"subject"],"rd.values.FoldEnrichment"])
NS.K27.fe <- cbind(NS.K27.overlap, NS.K27[NS.K27.overlap[,"subject"],"rd.values.FoldEnrichment"])
N.K4.fe <- cbind(N.K4.overlap, N.K4[N.K4.overlap[,"subject"],"rd.values.FoldEnrichment"])
N.K27.fe <- cbind(N.K27.overlap, N.K27[N.K27.overlap[,"subject"],"rd.values.FoldEnrichment"])

##get bin size based on virtual tss size and set all results tables to that size

res.NS.K4 <- res.NS.K27 <- res.N.K4 <- res.N.K27 <- integer(dim(v_tss)[1])

#sum the fold enrichment scores for each bin
nbins<- dim(v_tss)[1]
for(i in unique(NS.K4.fe[,1]) ){
 these <- NS.K4.fe[NS.K4.fe[,1]==i, 3]
 if(length(these)==0){
   res.NS.K4[i]<-0
 } else {
    res.NS.K4[i] <- sum(these)
  }
}

for(i in unique(NS.K27.fe[,1]) ){
 these <- NS.K27.fe[NS.K27.fe[,1]==i, 3]
 if(length(these)==0){
   res.NS.K27[i]<-0
 } else {
    res.NS.K27[i] <- sum(these)
  }
}

for(i in unique(N.K4.fe[,1]) ){
 these <- N.K4.fe[N.K4.fe[,1]==i, 3]
 if(length(these)==0){
   res.N.K4[i]<-0
 } else {
    res.N.K4[i] <- sum(these)
  }
}

for(i in unique(N.K27.fe[,1]) ){
 these <- N.K27.fe[N.K27.fe[,1]==i, 3]
 if(length(these)==0){
   res.N.K27[i]<-0
 } else {
    res.N.K27[i] <- sum(these)
  }
}



##plot it on a graph
v_tss.df <- as.data.frame(v_tss)

ind_neg <- which(v_tss.df[,"end"]<0)
ind_pos <- which(v_tss.df[,"end"]>0)

#take those indices and log them (removing negatives) them
v_tss.start.neg <- -(log10(-(v_tss.df[ind_neg,"start"])))
v_tss.start.pos <- log10(v_tss.df[ind_pos,"start"])

#then put them all back together
v_tss_start_log <- c(v_tss.start.neg,v_tss.start.pos)

postscript("dist_of_peaks_by_fold_enrichment.ps", horizontal = FALSE)
plot(v_tss.df[,"start"], 
     log2(res.NS.K4), 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="log(Total Fold Enrichment)",
     type="l"
)

lines(v_tss.df[,"start"], log2(res.NS.K27), col="red")
lines(v_tss.df[,"start"], log2(res.N.K4), col="green")
lines(v_tss.df[,"start"], log2(res.N.K27), col="magenta")

legend(75000, 120000, c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)
dev.off()

##########try plotting all peaks ordered by fold enrichment score
#take NS.K4.ranges file and add fold enrichments, order by FE and then plot using lines
NS.K4.ranges <- 


####################attempt to plot the same thing with larger p value cut offs to see if get more random peaks elsewhere
##cut off to only take the most significant peaks as loosing detail i think?!

NS.K4 <- NS.K4[order(NS.K4[,"rd.values.FoldEnrichment"],decreasing=TRUE),]
NS.K4 <- NS.K4[1:10,]
NS.K27 <- NS.K27[order(NS.K27[,"rd.values.FoldEnrichment"],decreasing=TRUE),]
NS.K27 <- NS.K27[1:10,]
N.K4 <- N.K4[order(N.K4[,"rd.values.FoldEnrichment"],decreasing=TRUE),]
N.K4 <- N.K4[1:10,]
N.K27 <- N.K27[order(N.K27[,"rd.values.FoldEnrichment"],decreasing=TRUE),]
N.K27 <- N.K27[1:10,]

library(IRanges)

minbp <--500
maxbp <-1000
bin.size <-5 
start <- seq(from=minbp, to=maxbp-bin.size, by=bin.size)
end   <- seq(from=minbp+bin.size-1, to=maxbp, by=bin.size)

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )

#calculate distances from peak start/end to tss
NS.K4.peakstart <- NS.K4[,"rd.start"] - NS.K4[,"tss.start"]
NS.K4.peakend <- NS.K4[,"rd.end"] - NS.K4[,"tss.start"]
NS.K4.ranges <- cbind(NS.K4.peakstart,NS.K4.peakend)
NS.K4.ir <- RangedData(ranges = IRanges(
                       start = NS.K4.ranges[,1],
                       end = NS.K4.ranges[,2],
                       )
                    )          

NS.K27.peakstart <- NS.K27[,"rd.start"] - NS.K27[,"tss.start"]
NS.K27.peakend <- NS.K27[,"rd.end"] - NS.K27[,"tss.start"]
NS.K27.ranges <- cbind(NS.K27.peakstart,NS.K27.peakend)
NS.K27.ir <- RangedData(ranges = IRanges(
                       start = NS.K27.ranges[,1],
                       end = NS.K27.ranges[,2],
                       )
                    )

N.K4.peakstart <- N.K4[,"rd.start"] - N.K4[,"tss.start"]
N.K4.peakend <- N.K4[,"rd.end"] - N.K4[,"tss.start"]
N.K4.ranges <- cbind(N.K4.peakstart,N.K4.peakend)
N.K4.ir <- RangedData(ranges = IRanges(
                       start = N.K4.ranges[,1],
                       end = N.K4.ranges[,2],
                       )
                    )

N.K27.peakstart <- N.K27[,"rd.start"] - N.K27[,"tss.start"]
N.K27.peakend <- N.K27[,"rd.end"] - N.K27[,"tss.start"]
N.K27.ranges <- cbind(N.K27.peakstart,N.K27.peakend)
N.K27.ir <- RangedData(ranges = IRanges(
                       start = N.K27.ranges[,1],
                       end = N.K27.ranges[,2],
                       )
                    )

#find overlaps of virtual tss 50bp regions with peaks
NS.K4.overlap <- findOverlaps(v_tss,
                              NS.K4.ir,
                              type = "within"
                              
                             )
NS.K4.overlap <- as.matrix(NS.K4.overlap[[1]])

NS.K27.overlap <- findOverlaps(v_tss, 
                              NS.K27.ir, 
                              type = "within", 
                             )
NS.K27.overlap <- as.matrix(NS.K27.overlap[[1]])

N.K4.overlap <- findOverlaps(v_tss, 
                              N.K4.ir, 
                              type = "within", 
                             )
N.K4.overlap <- as.matrix(N.K4.overlap[[1]])

N.K27.overlap <- findOverlaps(v_tss, 
                              N.K27.ir, 
                              type = "within", 
                             )
N.K27.overlap <- as.matrix(N.K27.overlap[[1]])

#count overlaps of virtual tss 50bp regions with peaks
NS.K4.count <- countOverlaps(v_tss, 
                              NS.K4.ir, 
                              type = "within", 
                             )
NS.K4.count <- as.data.frame(NS.K4.count)

NS.K27.count <- countOverlaps(v_tss, 
                              NS.K27.ir, 
                              type = "within", 
                             )
NS.K27.count <- as.data.frame(NS.K27.count)

N.K4.count <- countOverlaps(v_tss, 
                              N.K4.ir, 
                              type = "within", 
                             )
N.K4.count <- as.data.frame(N.K4.count)

N.K27.count <- countOverlaps(v_tss, 
                              N.K27.ir, 
                              type = "within", 
                             )
N.K27.count <- as.data.frame(N.K27.count)

v_tss.df <- as.data.frame(v_tss)

##############plot peak counts for each bin

postscript("peak_counts_v_tss_nonlog20-50k.ps", horizontal= FALSE)

plot(v_tss.df[,"start"], 
     (NS.K4.count[,"value"]),
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="Number of peaks",
     type="l"
    )
lines(v_tss.df[,"start"], (NS.K27.count[,"value"]), col="red")
lines(v_tss.df[,"start"], (N.K4.count[,"value"]), col="green")
lines(v_tss.df[,"start"], (N.K27.count[,"value"]), col="magenta")

legend(75000, 120000, c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)

dev.off()

##replace position in original (e.g. NS.K4) dataframe with "rd.values.FoldEnrichment"
NS.K4.fe <- cbind(NS.K4.overlap, NS.K4[NS.K4.overlap[,"subject"],"rd.values.FoldEnrichment"])
NS.K27.fe <- cbind(NS.K27.overlap, NS.K27[NS.K27.overlap[,"subject"],"rd.values.FoldEnrichment"])
N.K4.fe <- cbind(N.K4.overlap, N.K4[N.K4.overlap[,"subject"],"rd.values.FoldEnrichment"])
N.K27.fe <- cbind(N.K27.overlap, N.K27[N.K27.overlap[,"subject"],"rd.values.FoldEnrichment"])

##compact fold enrichments into one value for each position
bins <- unique(NS.K4.fe[,1])
res.NS.K4 <- integer(length(bins))

res.NS.K27 <- integer(length(bins))

res.N.K4 <- integer(length(bins))

res.N.K27 <- integer(length(bins))

for(i in bins){
 these <- NS.K4.fe[NS.K4.fe[,1]==i, 3]
 if (length(these) == 0) {
    res.NS.K4[i]=0 
 } else { 
    res.NS.K4[i] <- sum(these)/length(these)
 }
}

for(i in bins){
 these <- NS.K27.fe[NS.K27.fe[,1]==i, 3]
 res.NS.K27[i] <- sum(these)
}

for(i in bins){
 these <- N.K4.fe[N.K4.fe[,1]==i, 3]
 res.N.K4[i] <- sum(these)
}

for(i in bins){
 these <- N.K27.fe[N.K27.fe[,1]==i, 3]
 res.N.K27[i] <- sum(these)
}

##plot it on a graph somehow
v_tss.df <- as.data.frame(v_tss)

postscript("dist_of_peaks_by_fold_enrichment.ps", horizontal = FALSE)
plot(v_tss.df[,"start"], 
     res.NS.K4, 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="Cumulative Fold Enrichment",
     type="l"
)

lines(v_tss.df[,"start"], res.NS.K27, col="red")
lines(v_tss.df[,"start"], res.N.K4, col="green")
lines(v_tss.df[,"start"], res.N.K27, col="magenta")

legend(-10000, 60000, c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)







































dev.off()
























