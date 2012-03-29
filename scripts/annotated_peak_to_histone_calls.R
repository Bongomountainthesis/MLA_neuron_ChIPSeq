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





#######scripts below for merging histone data with expression data

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



########generate a smaller table with basic limma output and histone peak positions

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




########find bivalent genes and subset K4 and K27 genes
#bivalent.NS.inds <- intersect(which(hist.exp[,"NS_K4"]=="^chr.*"),(which(hist.exp[,"NS_K27"]=="^chr.*")))
####actually should be grep("^chr", hist.exp[,"NS_K4"]) and use value=T to get back values


##bivalent genes in NS and N
bivalent.NS.inds <- intersect(which(hist.exp[,"NS_K4"]!="<NA>"),which(hist.exp[,"NS_K27"]!="<NA>"))
bivalent.N.inds <- intersect(which(hist.exp[,"N_K4"]!="<NA>"),which(hist.exp[,"N_K27"]!="<NA>"))

bivalent.NS.genes <- hist.exp[bivalent.NS.inds,]
bivalent.N.genes <- hist.exp[bivalent.N.inds,]


##bivalent domains kept, lost or gained
bivalent.NS.kept <- which(bivalent.NS.genes[,"EnsemblID"] %in% bivalent.N.genes[,"EnsemblID"])
bivalent.NS.lost <- which(!(bivalent.NS.genes[,"EnsemblID"] %in% bivalent.N.genes[,"EnsemblID"]))
bivalent.N.gain <- which(!(bivalent.N.genes[,"EnsemblID"] %in% bivalent.NS.genes[,"EnsemblID"]))

bivalent.NS.kept.genes <- bivalent.NS.genes[bivalent.NS.kept,]
bivalent.NS.lost.genes <- bivalent.NS.genes[bivalent.NS.lost,]
bivalent.N.gain.genes <- bivalent.N.genes[bivalent.N.gain,]


##lost bivalent domains resolve to what?
bivalent.NS.lost.K4 <- which(bivalent.NS.lost.genes[,"N_K4"]!="<NA>")
bivalent.NS.lost.K27 <- which(bivalent.NS.lost.genes[,"N_K27"]!="<NA>")
bivalent.NS.lost.nothing <- intersect(which(is.na(bivalent.NS.lost.genes[,"N_K4"])), which(is.na(bivalent.NS.lost.genes[,"N_K27"])))

bivalent.NS.lost.K4.genes <- bivalent.NS.lost.genes[bivalent.NS.lost.K4,]
bivalent.NS.lost.K27.genes <- bivalent.NS.lost.genes[bivalent.NS.lost.K27,]
bivalent.NS.lost.nothing.genes <- bivalent.NS.lost.genes[bivalent.NS.lost.nothing,]


##write all to csv files

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

bivalent.NS.lost.K4.genes <- bivalent.NS.lost.K4.genes[order(bivalent.NS.lost.K4.genes[,"logFC"],decreasing=TRUE),]
write.csv(bivalent.NS.lost.K4.genes, "bivalent_NS_lost_to_K4_genes.csv")

bivalent.NS.lost.K27.genes <- bivalent.NS.lost.K27.genes[order(bivalent.NS.lost.K27.genes[,"logFC"],decreasing=TRUE),]
write.csv(bivalent.NS.lost.K27.genes, "bivalent_NS_lost_K27_genes.csv")

bivalent.NS.lost.nothing.genes <- bivalent.NS.lost.nothing.genes[order(bivalent.NS.lost.nothing.genes[,"logFC"], decreasing=TRUE),]
write.csv(bivalent.NS.lost.nothing.genes, "bivalent_NS_lost_nothing_genes.csv")
