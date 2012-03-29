#!/usr/bin/Rscript

#####################################################################################################
##script to remap all mohn data to mm9 (from mm8) and then link back to my data to allow comparison##
#####################################################################################################

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

chain.file <- "scripts/mm8ToMm9.over.chain"

qw <- function(...) {
  as.character(sys.call()[-1])
}


liftOver<-function(data, chain.file, ucsc.format=F, chr.col="chr", start.col="start",end.col="end"){

  #data should be a matrix or dataframe with cols for chr, start and end
  #TODO: Or a RangedData / IRange object, 
 
  this<-data.frame(chr=as.character(data[,chr.col]),start=data[,start.col],end=data[,end.col], stringsAsFactors=F)

  #Normal counting specifies ranges in 1-based, fully closed form.
  #UCSC specifies ranges in  0-based, half open
  if (!ucsc.format){
      this$start<-this$start-1
  }

  #use the chrstartend as an id so we can map back to the original data
  ids = paste(as.character(data[,chr.col]),data[,start.col],data[,end.col], sep=".")
  this <- cbind(this, ids)
  #If we have duplicate positions, remove them for the liftOver 
  this <- unique(this)
  #we also need to chuck out anything that doesn't have positional data.
  no.chr <- which(is.na(this[,1]))
  no.start <- which(is.na(this[,2]))
  no.end <- which(is.na(this[3]))

  inds <- unique(c(no.chr,no.start,no.end))

  if( length(inds)>0 ){ this <- this[-1*inds,] }
  
  
  ##all this stuff should be a .C() call but I don't have time to make it work just now.
  in.bed <- tempfile()
  write.table(this, file=in.bed, sep="\t", row.names=F, col.names=F, quote=F)

  out.bed <- tempfile()
  out.um.bed <- tempfile()
  lo.cmd <- paste("liftOver", in.bed, chain.file, out.bed, out.um.bed)
  system(lo.cmd)

  try(
    new.bed<-read.table(out.bed, sep="\t") 
      ,silent=T)
  try(
    new.um <- read.table(out.um.bed, sep="\t")
      ,silent=T)
  
  #throw away the files
  unlink(c(in.bed, out.bed, out.um.bed))

  if (!exists('new.bed')){stop("No successful mappings")}

  #use the ids as rownames chuck out the column
  #order in the same order as our original data
  #which should stick NAs in for anything that can't be mapped
  rownames(new.bed) <- new.bed[,4]
  new.bed <- new.bed[,-4]
  new.bed <- new.bed[ids,]

  if(!ucsc.format){
   #put the data back to 1-based
   new.bed[,2] <- new.bed[,2]+1
  }

  #replace the new positions in the original dataset
  data[,chr.col] <- new.bed[,1]
  data[,start.col] <- new.bed[,2]
  data[,end.col] <- new.bed[,3]

  #TODO: return some information about the data that won't map
  
  return(data)
  
}


##import data


library(ChIPpeakAnno)
library(biomaRt)

mohn <- read.csv("results/Mohn_2008_promoters.txt", sep="\t")

histone <- read.csv("results/histone_express_calls.csv")

##get chromosome start ends

position <- mohn[,"chromosomal_position_of_the_900bp_window_used_for_analysis_.MM8."]
position <- grep('^chr.*', position, value=T)

chr<-gsub('^(chr.*):.*','\\1', position)
start<-as.numeric(gsub('.*:(.*)-.*','\\1', position))
end<-as.numeric(gsub('.*-(.*)','\\1', position))

old.pos<-cbind(chr=chr, start=start, end=end)
new.pos<-liftOver(old.pos, chain.file=chain.file, ucsc.format=T)

old.pos.id <- paste(old.pos[,1],paste(old.pos[,2],old.pos[,3], sep="-"),sep=":")

new.pos <- cbind(new.pos, old.pos.id)

new.pos <- new.pos[which(!is.na(new.pos[,1])),]
new.pos[,"chr"] <- gsub('chr','',new.pos[,"chr"])

position.rd <- RangedData(ranges = IRanges(
                          start= as.numeric(new.pos[,"start"]),
                          end = as.numeric(new.pos[,"end"]),
                          ),
                          space = as.character(new.pos[,"chr"]),
                          values = new.pos[,"old.pos.id"]
                           )


#try straight from tss data

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
tss <- getAnnotation(ensmart, "TSS")


## problems with strandedness in ChIPPeakAnno? Do overlaps manually.
#convert tss data to a dataframe and then split into two to give strandedness
tss.df <- as.data.frame(tss)
tss.pos <- tss.df[which(tss.df[,"strand"]==1),]
tss.neg <- tss.df[which(tss.df[,"strand"]==-1),]
tss.pos <- RangedData(tss.pos)
tss.neg <- RangedData(tss.neg)


overlapping.tss.pos <- findOverlaps(position.rd,tss.pos)
overlapping.tss.neg <- findOverlaps(position.rd,tss.neg)

#turn everything into dataframe
position.rd.df <- as.data.frame(position.rd)
tss.pos.df <- as.data.frame(tss.pos)
tss.neg.df <- as.data.frame(tss.neg)

##process positive and negative separately and combine back together at the end somehow.

tss.olmatrix.pos <- as.matrix(overlapping.tss.pos)
tss.olmatrix.neg <- as.matrix(overlapping.tss.neg)

position.ind.pos <- tss.olmatrix.pos[,1]
position.ind.neg <- tss.olmatrix.neg[,1]
tss.ind.pos <- tss.olmatrix.pos[,2]
tss.ind.neg <- tss.olmatrix.neg[,2]

position.ol.pos <- position.rd.df[position.ind.pos,]
colnames(position.ol.pos) <- paste("rd", colnames(position.rd.df), sep=".")

position.ol.neg <- position.rd.df[position.ind.neg,]
colnames(position.ol.neg) <- paste("rd", colnames(position.rd.df), sep=".")

tss.pos.ol<-tss.pos.df[tss.ind.pos,]
colnames(tss.pos.ol) <- paste("tss", colnames(tss.pos.df), sep=".")

tss.neg.ol<-tss.neg.df[tss.ind.neg,]
colnames(tss.neg.ol) <- paste("tss", colnames(tss.neg.df), sep=".")

#combine rd.ol and tss data
res.pos <- cbind(position.ol.pos, tss.pos.ol)
inds.pos <- intersect(which(res.pos[,"tss.start"] > res.pos[,"rd.start"]), which(res.pos[,"rd.end"] > res.pos[,"tss.start"]))
res.pos<-res.pos[inds.pos,]

res.neg <- cbind(position.ol.neg, tss.neg.ol)
inds.neg <- intersect(which(res.neg[,"tss.end"] > res.neg[,"rd.start"]), which(res.neg[,"rd.end"] > res.neg[,"tss.end"]))
res.neg<-res.neg[inds.neg,]

#stick tables together
peaks <- rbind(res.pos, res.neg)

write.csv(peaks, "mohn_promoter_peaks_mm9.csv")

##remove all columns except rd.values (old pos name), ensembl ID and tss.description

keep.cols <- c("rd.values", "tss.names", "tss.description")

peaks <- peaks[,keep.cols]

#tidy up mohn data

keep.cols.mohn <- c(2, 3, 9, 24, 25, 26, 27, 28, 29)

mohn <- mohn[,keep.cols.mohn]

#merge data together

combine <- merge(mohn, 
                 peaks, 
                 by.x = "chromosomal_position_of_the_900bp_window_used_for_analysis_.MM8.",
                 by.y = "rd.values"
                )


#compare <- merge(histone, combine, by.x = "EnsemblID", by.y = "tss.names")
compare_all <- merge(histone, combine, by.x = "EnsemblID", by.y = "tss.names",all.x = TRUE,all.y = TRUE)

#order by K4+ and K27+ then remove duplicates - will get bivalent ones first if bivalent

compare <- compare_all[order(compare_all[,"K4_call_NP"], compare_all[,"K27_call_NP"],decreasing=TRUE),]

compare.nodups <- compare[!duplicated(compare[,"EnsemblID"]),]

##but still only 10615 genes that map between the two datasets - maybe try directly from the refseq id and ensembl ids again to check

write.csv(compare.nodups, "results/mine_vs_mohn.csv")

##work out which are bivalent in mine and mohn

mine.N.biv <- intersect(which(is.na(compare.nodups[,"N_K4"])),which(is.na(compare.nodupsres[,"N_K27"])))

mine.NS.biv <- intersect(which(res[,"NS_K4"]!="<NA>"),which(res[,"NS_K27"]!="<NA>"))

mohn.NS.biv <- intersect(which(res[,"K4_call_NP"]=="K4+"),which(res[,"K27_call_NP"]=="K27+"))

both.NS.biv <- intersect(intersect(which(res[,"NS_K4"]!="<NA>"),which(res[,"NS_K27"]!="<NA>")),intersect(which(res[,"K4_call_NP"]=="K4+"),which(res[,"K27_call_NP"]=="K27+")))





















##569 peaks have duplicate ensembl IDs - work out if any of them are bivalent



bivalent <- intersect(which(combine[,"K4_call_NP"]=="K4+"), which(combine[,"K27_call_NP"]=="K27+"))



biv <- c(1:length(combine[,1]))

combine <- cbind(combine, biv)

for(i in (combine[,"biv"]) ){
 combine[i,"biv"] <- "Mohn Bivalent Peak"
 }

bivalent <- combine[bivalent,c("tss.names", "biv")]

#take bivalent genes from mohn and merge back to my data

compare <- merge(histone, bivalent, by.x = "EnsemblID", by.y = "tss.names", all.x = TRUE)

d



