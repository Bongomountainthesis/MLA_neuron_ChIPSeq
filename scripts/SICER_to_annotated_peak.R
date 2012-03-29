#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

##########first load data in and make ranged data object from it

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
#filename = "MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups-W200-G200-islands-summary-FDR1E-3"

data <- read.csv(filename, header=F, sep="\t", comment.char = "#")

colnames(data) <- qw(Chr, Start, End, nTags_ChIP, nTags_Cnt, PValue, FoldEnrichment, FDR)

##calculate neg10 log pvalues
log <- -10*(log10(data[,"PValue"]))

##calculate length of peak
length <- data[,"End"] - data[,"Start"]

##add to dataframe
data_out <- data.frame(data,log, length)
colnames(data_out) <- qw(Chr, Start, End, nTags_ChIP, nTags_Cnt, PValue, FoldEnrichment, FDR, neg10log10pVal, Length)

###remove peaks with crap FDRs - cut at <1e-20 for the time being

data_out <- data_out[data_out[,"FDR"]<3E-8,]

values <-  qw(Length, nTags_ChIP, nTags_Cnt, neg10log10pVal, FDR, FoldEnrichment)

#sometimes they have FDR
if(ncol(data)==9){
  colnames(data)[9] <- "FDR"
  values<-c(values,"FDR")
}

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
data_out[,"Name"]<-paste(paste(data_out[,"Chr"],data_out[,"Start"],sep=":"),data_out[,"End"], sep="-")


#note - this fails if chr position is undefined, which can happen
#if we've mapped over from somewhere else, so:
ids <- which(is.na(data[,"Chr"]) | is.na(data[,"Start"]) | is.na(data[,"End"]) )
if(length(ids)>0){
  data <- data[-ids,]
}

rd <- RangedData(ranges = IRanges(
                   start= data_out$Start,
                   end = data_out$End,
                   names = as.character(data_out$Name),
                   ),
                 space = as.character(data_out$Chr),
                 values = data_out[,values]
                 )

#outfile=sub('.','.RangedData.RData',filename)

#And save the result
save(rd, file="./SICER_peaks.RangedData.RData")






##########Now annotate peaks with overlapping tss start

library(ChIPpeakAnno)
library(biomaRt)
#ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


#load mouse transcripts
data(TSS.mouse.NCBIM37)

#The chromosomes names used in teh NCBIM37 dataset are 1..Y, MT and
#various NT_***** non-standard chrs. 

                                        #remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}

#chippeakanno discards anything in values, so hang on to them
#vals <- as.data.frame(rd)
#vals <- vals[,grep("values.", colnames(vals), value=T)]
#colnames(vals) <- gsub('values.','',colnames(vals))


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions,
# not individual transcripts.

#get the most recent annotation data from ensembl
#tss <- getAnnotation(ensmart, "TSS")
#save(tss, file="tss.RData")
tss <- get(load("tss.RData"))

#convert tss data to a dataframe and then split into two to give strandedness
tss.df <- as.data.frame(tss)
tss.pos <- tss.df[which(tss.df[,"strand"]==1),]
tss.neg <- tss.df[which(tss.df[,"strand"]==-1),]
tss.pos <- RangedData(tss.pos)
tss.neg <- RangedData(tss.neg)


#and later we'll want to know if it overlaps a coding region?
#exons <- getAnnotation(ensmart, "Exon")
#save(exons, file="exons.RData")
#exons <- get(load("exons.RData"))

#rd.df <- data.frame(rd)


#### try findOverlaps function instead as ChIPPeakAnno is a big pile of arse
#find overlaps on the positive strand

overlapping.tss.pos <- findOverlaps(rd,tss.pos)
overlapping.tss.neg <- findOverlaps(rd,tss.neg)

#turn everything into dataframe
rd.df <- as.data.frame(rd)
tss.pos.df <- as.data.frame(tss.pos)
tss.neg.df <- as.data.frame(tss.neg)

##process positive and negative separately and combine back together at the end somehow.

tss.olmatrix.pos <- as.matrix(overlapping.tss.pos)
tss.olmatrix.neg <- as.matrix(overlapping.tss.neg)

rd.ind.pos <- tss.olmatrix.pos[,1]
rd.ind.neg <- tss.olmatrix.neg[,1]
tss.ind.pos <- tss.olmatrix.pos[,2]
tss.ind.neg <- tss.olmatrix.neg[,2]

rd.ol.pos <- rd.df[rd.ind.pos,]
colnames(rd.ol.pos) <- paste("rd", colnames(rd.df), sep=".")

rd.ol.neg <- rd.df[rd.ind.neg,]
colnames(rd.ol.neg) <- paste("rd", colnames(rd.df), sep=".")

tss.pos.ol<-tss.pos.df[tss.ind.pos,]
colnames(tss.pos.ol) <- paste("tss", colnames(tss.pos.df), sep=".")

tss.neg.ol<-tss.neg.df[tss.ind.neg,]
colnames(tss.neg.ol) <- paste("tss", colnames(tss.neg.df), sep=".")

#combine rd.ol and tss data
res.pos <- cbind(rd.ol.pos, tss.pos.ol)
inds.pos <- intersect(which(res.pos[,"tss.start"] > res.pos[,"rd.start"]), which(res.pos[,"rd.end"] > res.pos[,"tss.start"]))
res.pos<-res.pos[inds.pos,]

res.neg <- cbind(rd.ol.neg, tss.neg.ol)
inds.neg <- intersect(which(res.neg[,"tss.end"] > res.neg[,"rd.start"]), which(res.neg[,"rd.end"] > res.neg[,"tss.end"]))
res.neg<-res.neg[inds.neg,]

#stick tables together
peaks <- rbind(res.pos, res.neg)

write.csv(peaks, "all_overlapping_peaks.csv")



