
###########take raw expression data from BSData and link back to histone calls data and then compare actual gene expression levels to histone calls

library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)

BSData <- get(load("data/BSData.quantile.RData"))
E <- exprs(BSData)

hist.exp <- read.csv("results/histone_express_calls.csv")

####add Ensembl IDs to E as something wierd is happening....

ids = rownames(E)
symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, illuminaMousev2BeadIDENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1) 
length(crosshyb)
ensembl[crosshyb] <- NA 
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

##########make averages for NS and N - use apply function - before just stuck in excel as for loops didnt work
NS.avg <-integer(dim(E)[1])

ns.cols <- grep("NS", colnames(E))
n.cols  <- grep("N[^S]", colnames(E)) 

aves <- apply(E, 1, function(x){  
    NS.ave  <- sum(x[ns.cols])/length(ns.cols)
    N.ave   <- sum(x[n.cols])/length(n.cols)
    return(c(NS.ave, N.ave))
 }
)

aves <- t(aves)
colnames(aves)<-c("NS.ave", "N.ave")

hist.exp <- cbind(hist.exp, aves[as.character(hist.exp[,"ID"]),])
hist.exp <- hist.exp[,-1]


####sod this....for loops hate me - only works for some of the columns but not all. Make csv file and stick in excel for the time being


###add back to E to give ensembl ids and raw expression values

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   E[,1:8],
                   stringsAsFactors=F
              )

#write.csv(anno, "raw_probe_expression_annotated.csv")

#anno <- read.csv("raw_probe_expression_annotated.csv")

#####remove probes in anno that are not in hist.exp (ie the duplicated probes)

inds <- which(anno[,"ID"] %in% hist.exp[,"ID"])

anno.nodups <- anno[inds,]

##tidy it up abit before merging...

colrm.anno.nodups <- c("X", "ID")
keep.cols <- !(colnames(anno.nodups) %in% colrm.anno.nodups)
anno.nodups <- anno.nodups[,keep.cols]

hist <- merge(hist.exp, anno.nodups, by.x = "EnsemblID", by.y = "EnsemblID")

write.csv(hist, "hist_exp_plus_rawexpvalues.csv")

















