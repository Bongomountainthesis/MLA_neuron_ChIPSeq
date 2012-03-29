#!/usr/local/bin/Rscript


###load limma results file and nearest or overlapping feature data and merge them together

limma.filename <- "/space/matt/mla_neurons_lumixpn/results/limma_results.csv"
NS.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K4me3/overlapping_peaks.csv"
NS.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K27me3/overlapping_peaks.csv"
N.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K4me3/overlapping_peaks.csv"
N.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K27me3/overlapping_peaks.csv"

limma <- read.csv(limma.filename)
NS.K4 <- read.csv(NS.K4.filename)
NS.K27 <- read.csv(NS.K27.filename)
N.K4 <- read.csv(N.K4.filename)
N.K27 <- read.csv(N.K27.filename)


#calculate distance from TSS for each peak (middle) - this below ignores strandedness DO NOT USE

NS.K4.dist <- (((NS.K4[,"rd.end"] - NS.K4[,"rd.start"])/2) + NS.K4[,"rd.start"]) - NS.K4[,"tss.start"]
NS.K27.dist <- (((NS.K27[,"rd.end"] - NS.K27[,"rd.start"])/2) + NS.K27[,"rd.start"]) - NS.K27[,"tss.start"]
N.K4.dist <- (((N.K4[,"rd.end"] - N.K4[,"rd.start"])/2) + N.K4[,"rd.start"]) - N.K4[,"tss.start"]
N.K27.dist <- (((N.K27[,"rd.end"] - N.K27[,"rd.start"])/2) + N.K27[,"rd.start"]) - N.K27[,"tss.start"]

NS.K4 <- cbind(NS.K4, NS.K4.dist)
NS.K27 <- cbind(NS.K27, NS.K27.dist)
N.K4 <- cbind(N.K4, N.K4.dist)
N.K27 <- cbind(N.K27, N.K27.dist)

##remove duplicate probe-gene matches - keep highest FC difference
limma.o <- limma[order(abs(limma[,"logFC"]), decreasing=TRUE),]
limma.do <- limma.o[!duplicated(limma.o[,"EnsemblID"]),]

postscript(file="results/sigvsdist.ps", horizontal=FALSE)
plot(NS.K4[,"NS.K4.dist"],
     NS.K4[,"rd.values.neg10log10pVal"], 
     xlim=c(-10000,10000), 
     pch=".", 
     main="Distribution of Peaks (middle) around TSS", 
     ylab="Fold enrichment", 
     xlab="Distance to TSS (bps)"
    )
    points(NS.K27[,"NS.K27.dist"], NS.K27[,"rd.values.neg10log10pVal"], pch=".", col = "red")
    points(N.K4[,"N.K4.dist"], N.K4[,"rd.values.neg10log10pVal"], pch=".", col = "blue")
    points(N.K27[,"N.K27.dist"], N.K27[,"rd.values.neg10log10pVal"], pch=".", col = "darkgreen")

    #abline(v=4000, col="red")
    #abline(v=-4000, col="red")
    #abline(h=100, col="red")
dev.off()

library(IRanges)

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
library(IRanges)

minbp <--50000
maxbp <-50000
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



##plot it on a graph - add details of longest sets of peaks (ie NS.K27) to bins to give all peaks across range with the same lengths of x for plot
v_tss.df <- as.data.frame(v_tss)

ind_neg <- which(v_tss.df[,"end"]<0)
ind_pos <- which(v_tss.df[,"end"]>0)
#take those indices and postive them
v_tss.start.neg <- -(log2(-(v_tss.df[ind_neg,"start"])))

#then put them all back together
v_tss.start.pos <- log2(v_tss.df[ind_pos,"start"])

v_tss_start_log <- c(v_tss.start.neg,v_tss.start.pos)

postscript("dist_of_peaks_by_fold_enrichment.ps", horizontal = FALSE)
plot(v_tss_start_log, 
     log2(res.NS.K4), 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="log2(Fold Enrichment)",
     type="l"
)

lines(v_tss_start_log, log2(res.NS.K27), col="red")
lines(v_tss_start_log, log2(res.N.K4), col="green")
lines(v_tss_start_log, log2(res.N.K27), col="magenta")

legend(-50000, 17, c("NS-K4","NS-K27", "N-K4", "N-K27"), 
   col=c("blue","red","green","magenta"), pch=21:22, lty=1:2)
dev.off()

##########try plotting all peaks ordered by fold enrichment score
#take NS.K4.ranges file and add fold enrichments, order by FE and then plot using lines
NS.K4.ranges <- 


####################attempt to plot the same thing with larger p value cut offs to see if get more random peaks elsewhere
library(IRanges)

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


##merge histone calls to expression data

#add tag to each column name
colnames(NS.K4) <- paste("NS_K4", colnames(NS.K4),sep="_")
colnames(NS.K27) <- paste("NS_K27", colnames(NS.K27),sep="_")
colnames(N.K4) <- paste("N_K4", colnames(N.K4),sep="_")
colnames(N.K27) <- paste("N_K27", colnames(N.K27),sep="_")

##merge columns by ensembl id
dat.NS.K4 <- merge(limma.do, NS.K4, by.x="EnsemblID", by.y="NS_K4_tss.names", all.x=T)
dat.NS.K27 <- merge(dat.NS.K4, NS.K27, by.x="EnsemblID", by.y="NS_K27_tss.names", all.x=T)
dat.N.K4 <- merge(dat.NS.K27, N.K4, by.x="EnsemblID", by.y="N_K4_tss.names", all.x=T)
dat.N.K27 <- merge(dat.N.K4, N.K27, by.x="EnsemblID", by.y="N_K27_tss.names", all.x=T)

dat.all <- dat.N.K27

write.csv(dat.all, "/space/matt/MLA_neuron_ChIPSeq/results/expression_histones_all.csv")



##generate table with Ensembl ID, gene symbol and description, logFC, adj.PVal and whether it 'has' K4 and K27 in NS and Neuron
hist <- data.frame(dat.all[,"NS_K4_rd.names"],dat.all[,"NS_K27_rd.names"],dat.all[,"N_K4_rd.names"],dat.all[,"N_K27_rd.names"])
colnames(hist) <- c("NS_K4", "NS_K27", "N_K4", "N_K27")

##order limma.do by EnsemblID so that making a dataframe puts the columns in the same order
limma.doo <- limma.do[order(limma.do[,"EnsemblID"], decreasing=FALSE),]

hist.exp <- data.frame(limma.doo, hist)
rm.cols<-c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "AveExpr", "t", "P.Value", "B")
keep.cols<-!(colnames(hist.exp) %in% rm.cols)
hist.exp <- hist.exp[,keep.cols]

write.csv(hist.exp, "/space/matt/MLA_neuron_ChIPSeq/results/histone_express_calls.csv")

###find bivalent genes and subset K4 and K27 genes
#bivalent.NS.inds <- intersect(which(hist.exp[,"NS_K4"]=="^chr.*"),(which(hist.exp[,"NS_K27"]=="^chr.*")))
####actually should be grep("^chr", hist.exp[,"NS_K4"]) and use value=T to get back values



bivalent.NS.inds <- intersect(which(hist.exp[,"NS_K4"]!="<NA>"),which(hist.exp[,"NS_K27"]!="<NA>"))
bivalent.N.inds <- intersect(which(hist.exp[,"N_K4"]!="<NA>"),which(hist.exp[,"N_K27"]!="<NA>"))

bivalent.NS.genes <- hist.exp[bivalent.NS.inds,]
bivalent.N.genes <- hist.exp[bivalent.N.inds,]

bivalent.NS.kept <- which(bivalent.NS.genes[,"EnsemblID"] %in% bivalent.N.genes[,"EnsemblID"])
bivalent.NS.lost <- which(!(bivalent.NS.genes[,"EnsemblID"] %in% bivalent.N.genes[,"EnsemblID"]))
bivalent.N.gain <- which(!(bivalent.N.genes[,"EnsemblID"] %in% bivalent.NS.genes[,"EnsemblID"]))

bivalent.NS.kept.genes <- bivalent.NS.genes[bivalent.NS.kept,]
bivalent.NS.lost.genes <- bivalent.NS.genes[bivalent.NS.lost,]
bivalent.N.gain.genes <- bivalent.N.genes[bivalent.N.gain,]

bivalent.NS.genes <- bivalent.NS.genes[order(bivalent.NS.genes[,"logFC"], decreasing=TRUE),]
write.csv(bivalent.NS.genes,"bivalent_NS_genes.csv")

bivalent.N.genes <- bivalent.N.genes[order(bivalent.N.genes[,"logFC"], decreasing=TRUE),]
write.csv(bivalent.N.genes,"bivalent_N_genes.csv")

bivalent.NS.kept.genes <- bivalent.NS.kept.genes[order(bivalent.NS.kept.genes[,"logFC"], decreasing=TRUE),]
write.csv(bivalent.NS.kept.genes, "bivalent_NS_kept_genes.csv")

bivalent.NS.lost.genes <- bivalent.NS.lost.genes[order(bivalent.NS.lost.genes[,"logFC"],decreasing=TRUE),]
write.csv(bivalent.NS.lost.genes,"bivalent_NS_lost_genes.csv")

bivalent.N.gain.genes <- bivalent.N.gain.genes[order(bivalent.N.gain.genes[,"logFC"],decreasing=TRUE),]
write.csv(bivalent.N.gain.genes,"bivalent_N_gain_genes.csv")




#####################################turn histone calls into + or -

hist.exp <- factor(as.character(gsub("TRUE", "chr.*", hist.exp[,7])))

hist.exp <- data.frame(limma.do, hist)

#remove all the spare columns

rm.cols<-c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "AveExpr", "t", "P.Value", "B")
keep.cols<-!(colnames(hist.exp) %in% rm.cols)
hist.exp <- hist.exp[,keep.cols]

#turn TRUE into + and NA into -

NSK4 <- factor(as.character(gsub("TRUE","+",hist.exp[,"NS_K4"])))
NSK27 <- factor(as.character(gsub("TRUE","+",hist.exp[,"NS_K27"])))
NK4 <- factor(as.character(gsub("TRUE","+",hist.exp[,"N_K4"])))
NK27 <- factor(as.character(gsub("TRUE","+",hist.exp[,"N_K27"])))

hist.calls <- data.frame(NSK4, NSK27, NK4, NK27)
#hist.calls <- na.omit(hist.calls)


hist.calls.exp <- data.frame(limma.do,hist.calls)
keep.cols<-!(colnames(hist.calls.exp) %in% rm.cols)
hist.exp <- hist.calls.exp[,keep.cols]

write.csv(hist.exp,"/space/matt/MLA_neuron_ChIPSeq/results/histone_expression_calls.csv")





















