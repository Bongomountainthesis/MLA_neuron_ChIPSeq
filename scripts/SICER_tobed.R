#!/usr/local/bin/Rscript

dat <- read.csv("nearest_tss_all.csv")

dat.p <- dat[dat[,"FDR"]<1E-20,]

#colnms <- c("Chr", "Start", "End", "Name", "Score")

bed <- dat.p[,c("space", "start", "end", "names", "neg10log10pVal")]
bed[,1] <- paste("chr", bed[,1], sep="")

write.table(bed, file="peaks.bed", col.names=F, row.names=F, sep="\t", quote=T)

