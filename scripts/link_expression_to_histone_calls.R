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

#calculate distance from TSS for each peak (middle)
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
#split virtual tss into 100bp segments across 2kb and then use findoverlaps to give a density plot thing
v_tss <- read.csv("/space/matt/MLA_neuron_ChIPSeq/results/virtual_tss_ps.csv")
v_tss <- RangedData(ranges = IRanges(
                   start= v_tss$Start,
                   end = v_tss$End,
                   names = as.character(v_tss$Name),
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

#count overlaps of virtual tss 50bp regions with peaks
NS.K4.overlap <- countOverlaps(v_tss, 
                              NS.K4.ir, 
                              type = "within", 
                             )
NS.K4.overlap <- as.data.frame(NS.K4.overlap)
#NS.K4.overlap <- cbind(c(1:160), NS.K4.overlap)

NS.K27.overlap <- countOverlaps(v_tss, 
                              NS.K27.ir, 
                              type = "within", 
                             )
NS.K27.overlap <- as.data.frame(NS.K27.overlap)

N.K4.overlap <- countOverlaps(v_tss, 
                              N.K4.ir, 
                              type = "within", 
                             )
N.K4.overlap <- as.data.frame(N.K4.overlap)

N.K27.overlap <- countOverlaps(v_tss, 
                              N.K27.ir, 
                              type = "within", 
                             )
N.K27.overlap <- as.data.frame(N.K27.overlap)

v_tss.df <- as.data.frame(v_tss)

postscript("dist_of_peaks.ps", horizontal = FALSE)
plot(v_tss.df[,"start"], 
     NS.K4.overlap[,"value"], 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="Number of peaks"
)

points(v_tss.df[,"start"], NS.K27.overlap[,"value"], col="red")
points(v_tss.df[,"start"], N.K4.overlap[,"value"], col="blue", pch="x")
points(v_tss.df[,"start"], N.K27.overlap[,"value"], col="red", pch="x")

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





















