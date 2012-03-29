#!/usr/local/bin/Rscript


###load limma results file and nearest or overlapping feature data and merge them together

limma.filename <- "/space/matt/mla_neurons_lumixpn/results/limma_results.csv"
NS.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K4me3/nearest_tss_NS_H3K4me3.csv"
NS.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K27me3/nearest_tss_NS_H3K27me3.csv"
N.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K4me3/nearest_tss_N_H3K4me3.csv"
N.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K27me3/nearest_tss_N_H3K27me3.csv"

limma <- read.csv(limma.filename)
NS.K4 <- read.csv(NS.K4.filename)
NS.K27 <- read.csv(NS.K27.filename)
N.K4 <- read.csv(N.K4.filename)
N.K27 <- read.csv(N.K27.filename)

##remove duplicate probe-gene matches - keep highest FC difference
limma.o <- limma[order(abs(limma[,"logFC"]), decreasing=TRUE),]
limma.do <- limma.o[!duplicated(limma.o[,"EnsemblID"]),]

##calculate which genes have a 'real' peak
#for K4 cut off at 1kb and neg10logpVal at 250
#for K27 cut off at 3kb
postscript(file="results/NS_K27_sigvsdist.ps", horizontal=FALSE)
plot(NS.K27[,"distancetoFeature"],
     NS.K27[,"neg10log10pVal"], 
     xlim=c(-10000,10000), 
     pch=".", 
     main="Distribution of Peaks around TSS", 
     ylab="-10*(log10PValue)", 
     xlab="Distance to TSS (bps)"
    )
    abline(v=4000, col="red")
    abline(v=-4000, col="red")
    #abline(h=100, col="red")
dev.off()

NS.K4.dist <- NS.K4[abs(NS.K4[,"distancetoFeature"])<=1500,]
NS.K4.sig <- NS.K4.dist[NS.K4.dist[,"neg10log10pVal"]>=0,]

NS.K27.dist <- NS.K27[abs(NS.K27[,"distancetoFeature"])<=4000,]
NS.K27.sig <- NS.K27.dist[NS.K27.dist[,"neg10log10pVal"]>=0,]

N.K4.dist <- N.K4[abs(N.K4[,"distancetoFeature"])<=1500,]
N.K4.sig <- N.K4.dist[N.K4.dist[,"neg10log10pVal"]>=0,]

N.K27.dist <- N.K27[abs(N.K27[,"distancetoFeature"])<=4000,]
N.K27.sig <- N.K27.dist[N.K27.dist[,"neg10log10pVal"]>=0,]

# for each gene, keep the info only for the nearest peak
nms <- NS.K4.sig[,"feature"]
dis <- abs(NS.K4.sig[,"distancetoFeature"])
nearest.NS_K4 <- NS.K4.sig[order(nms, dis, decreasing=TRUE),]
nearest.NS_K4 <-nearest.NS_K4[!duplicated(nearest.NS_K4[,"feature"]),]

nms <- NS.K27.sig[,"feature"]
dis <- abs(NS.K27.sig[,"distancetoFeature"])
nearest.NS_K27 <- NS.K27.sig[order(nms, dis, decreasing=TRUE),]
nearest.NS_K27 <-nearest.NS_K27[!duplicated(nearest.NS_K27[,"feature"]),]

nms <- N.K4.sig[,"feature"]
dis <- abs(N.K4.sig[,"distancetoFeature"])
nearest.N_K4 <- N.K4.sig[order(nms, dis, decreasing=TRUE),]
nearest.N_K4 <-nearest.N_K4[!duplicated(nearest.N_K4[,"feature"]),]

nms <- N.K27.sig[,"feature"]
dis <- abs(N.K27.sig[,"distancetoFeature"])
nearest.N_K27 <- N.K27.sig[order(nms, dis, decreasing=TRUE),]
nearest.N_K27 <-nearest.N_K27[!duplicated(nearest.N_K27[,"feature"]),]


#add tag to each column name
colnames(nearest.NS_K4) <- paste("NS_K4", colnames(nearest.NS_K4),sep="_")
colnames(nearest.NS_K27) <- paste("NS_K27", colnames(nearest.NS_K27),sep="_")
colnames(nearest.N_K4) <- paste("N_K4", colnames(nearest.N_K4),sep="_")
colnames(nearest.N_K27) <- paste("N_K27", colnames(nearest.N_K27),sep="_")

##merge columns by ensembl id
dat.NS.K4 <- merge(limma.do, nearest.NS_K4, by.x="EnsemblID", by.y="NS_K4_feature", all.x=T)
dat.NS.K27 <- merge(dat.NS.K4, nearest.NS_K27, by.x="EnsemblID", by.y="NS_K27_feature", all.x=T)
dat.N.K4 <- merge(dat.NS.K27, nearest.N_K4, by.x="EnsemblID", by.y="N_K4_feature", all.x=T)
dat.N.K27 <- merge(dat.N.K4, nearest.N_K27, by.x="EnsemblID", by.y="N_K27_feature", all.x=T)

dat.all <- dat.N.K27

write.csv(dat.all, "/space/matt/MLA_neuron_ChIPSeq/results/expression_histones_all.csv")



##generate table with Ensembl ID, gene symbol and description, logFC, adj.PVal and whether it 'has' K4 and K27 in NS and Neuron
hist <- data.frame(dat.all[,"NS_K4_names"],dat.all[,"NS_K27_names"],dat.all[,"N_K4_names"],dat.all[,"N_K27_names"])
colnames(hist) <- c("NS_K4", "NS_K27", "N_K4", "N_K27")

##order limma.do by EnsemblID so that making a dataframe puts the columns in the same order
limma.doo <- limma.do[order(limma.do[,"EnsemblID"], decreasing=FALSE),]

hist.exp <- data.frame(limma.doo, hist)
rm.cols<-c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "AveExpr", "t", "P.Value", "B")
keep.cols<-!(colnames(hist.exp) %in% rm.cols)
hist.exp <- hist.exp[,keep.cols]

write.csv(hist.exp, "/space/matt/MLA_neuron_ChIPSeq/results/histone_express_calls.csv")




####################not sure this bit is doing what i think its supposed to be?!

hist <- cbind(dat.NS.K4[,"NS_K4_names"]!="<NA>",dat.NS.K27[,"NS_K27_names"]!="<NA>",dat.N.K4[,"N_K4_names"]!="<NA>",dat.N.K27[,"N_K27_names"]!="<NA>")
colnames(hist) <- c("NS_K4", "NS_K27", "N_K4", "N_K27")
hist <- data.frame(hist)



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





















