
options(stringAsFactors=FALSE)

### load in data (will take awhile)

#k4 <- read.csv(file = "MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups.bed",header = FALSE, sep = "\t")

#colnames(k4) <- c("Chr", "Start", "End", "Stuff", "Something", "Strand")

#save(k4, file = "k4_reads.RData")

##just for testing take 100 lines
#k4_full <- k4

#k4 <- k4[1:100,]

#library(IRanges)

### make ranged data object

## make peak name

#name <- paste(paste(k4[,"Chr"],k4[,"Start"],sep=":"),k4[,"End"], sep="-")

#data <- cbind(k4[,"Chr"],k4[,"Start"],k4[,"End"]

#rd <- RangedData(ranges = IRanges(
#			start = k4$Start,
#			end = k4$End,
#			names = as.character(k4$Name),
#			),
#			space = as.character(k4$Chr)
#		)

## annotate to nearest TSS

#library(ChIPpeakAnno)
#library(biomaRt)
#ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

#get the most recent annotation data from ensembl
#tss <- getAnnotation(ensmart, "TSS")
#save(tss, file=paste(dirname(filename),"/tss.RData",sep=""))

#peak_to_gene <- annotatePeakInBatch(rd,
#                                    AnnotationData=tss,
#                                   PeakLocForDistance = "middle",    # from the middle of the peak
#                                    FeatureLocForDistance = "TSS",  # to the TSS of the gene
#                                    output = "both",
#                                    multiple=TRUE
#                                    )



################################### 

##try again using bam reads that are in peaks

options(stringsAsFactors = FALSE);

library(biomaRt)
library(IRanges)
library(ShortRead)
library(snow)
library(baySeq)
library(DESeq)
library(rtracklayer)
library(Rsamtools)
library(ChIPpeakAnno)


###need to use ScanBamParam to pull "pos" from bam file
bamFile <- "MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups.bam"

### need to make rd object from bed file

#bed <- read.csv(file="SICER_MLA_H3K4me3/MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups-W200-G200-islands-summary-FDR1E-3",sep="\t",header = FALSE)

# try adding already selecting peaks based on overlapping of TSS
bed <- read.csv(file="all_overlapping_peaks.csv")
bed <- bed[,c(2,3,4,8,9,10,11)]

##take top 100 peaks

colnames(bed) <- c("Chr","Start","End","nTags_ChIP", "nTags_Cnt", "neg10log10pVal","FDR")

##need to add "chr" to columns...

bed[,"Chr"] <- paste("chr",bed[,"Chr"],sep = "")

##need to remove nonunique reads - not entirely sure why they're still in there??

bed.o <- bed[order(bed[,"neg10log10pVal"],decreasing=TRUE),]
#bed.top <- bed.o[!duplicated(bed.o[,"Start"]),]
bed.top <- bed.o[900,]

name <- paste(paste(bed.top[,"Chr"],bed.top[,"Start"],sep=":"),bed.top[,"End"], sep="-")

data <- data.frame(name,bed.top[,"Chr"],bed.top[,"Start"],bed.top[,"End"])

colnames(data) <- c("Name","Chr","Start","End")

rd <- RangedData(ranges = IRanges(
			start = data$Start,
			end = data$End,
			names = as.character(data$Name),
			),
			space = as.character(data$Chr)
		)

what <- c("rname","pos") 
param <- ScanBamParam(which=rd, what=what)
bam <- scanBam(bamFile, param=param)
bam.reads <- sapply(bam, function(x){paste(paste(x$rname,x$pos,sep=":"),x$pos+35,sep="-")})

####### now annotate them with nearest TSS start

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")

##add NAs to extra lists
#maxLen <- max(sapply(bam.reads, length))
# create a new list with elements padded out with NAs
#new.bam.reads <- lapply(bam.reads, function(.ele){
#    c(.ele, rep(NA, maxLen))[1:maxLen]
#  })

#bam.matrix <- do.call(rbind, new.bam.reads)

#bam.cat <- apply(bam.matrix, 2, rbind)


bam.unlist <- unlist(bam.reads)

###parse back out to make rd object

#take only unique reads

bam.res <- unique(bam.unlist)

res <- data.frame(integer(length(bam.res)-1),integer(length(bam.res)-1),integer(length(bam.res)-1),integer(length(bam.res)-1))

colnames(res) <- c("Name","Chr","Start","End")

##for some reason puts an extra row on
#res <- res[1:9360,]

for(i in 1:length(bam.res)){
		res[i,"Name"] <- bam.res[i]
		res[i,"Chr"] <- substr(bam.res[i],4,which(strsplit(bam.res[i],'')[[1]]==':')-1)
		res[i,"Start"] <- substr(gsub('.*\\:','',bam.res[i]),0, which(strsplit(gsub('.*\\:','',bam.res[i]),'')[[1]]=='-')-1)
		res[i,"End"] <- gsub('.*\\-','',bam.res[i])
		}

rd <- RangedData(ranges = IRanges(
			start = as.numeric(res$Start),
			end = as.numeric(res$End),
			names = as.character(res$Name),
			),
			space = as.character(res$Chr)
		)


       
peak_to_gene <- annotatePeakInBatch(rd,
                                    AnnotationData=tss,
                                    PeakLocForDistance = "start",    # from the middle of the peak
                                    FeatureLocForDistance = "TSS",    # to the TSS of the gene
                                    output = "both",
                                    multiple=TRUE
                                    )

##take out of rd

peak_to_gene.df <- as.data.frame(peak_to_gene)

##convert read start/end into versus TSS start

#first separate pos/neg strands

read.pos <- cbind(peak_to_gene.df[,"distancetoFeature"], (peak_to_gene.df[,"distancetoFeature"]+35))
colnames(read.pos) <- c("Start","End")

reads <- RangedData(ranges = IRanges(
			start = read.pos[,"Start"],
			end = read.pos[,"End"]
			))

###need to take peaks only over tss - ie go back to peak list and pull out those bam files...

#count overlaps of virtual tss 50bp regions with peaks

minbp <--5000
maxbp <-10000
bin.size <-10 
start <- seq(from=minbp, to=maxbp-bin.size, by=bin.size)
end   <- seq(from=minbp+bin.size-1, to=maxbp, by=bin.size)

v_tss <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   )
                 )

overlap <- findOverlaps(v_tss,
                              reads,
                              type = "within"
                              
                             )
read.overlap <- as.matrix(overlap[[1]])

count <- rle(read.overlap[,1])

query <- unique(read.overlap[,"query"])

res <- cbind(query, count$lengths)

tss.res <- data.frame(integer(length(start)),integer(length(start)))
tss.res[,1] <- 1:length(start)
tss.res[,2] <- start


colnames(tss.res) <- c("Bin","Position")
 
res.table <- merge(tss.res, res, by.x = "Bin", by.y = "query", all.x = TRUE)

res.table[is.na(res.table)] <- 0 

colnames(res.table) <- c("Bin","Position","Count")

postscript("dist_of_peaks_by_fold_enrichment.ps", horizontal = FALSE)
plot(res.table[,"Position"], 
     log2(res.table[,"Count"]), 
     col="blue", 
     main="Distribution of peaks around a virtual TSS", 
     xlab= "Distance to TSS", 
     ylab="Number of Reads)",
     type="l"
)

dev.off()
