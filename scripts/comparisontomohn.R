options(stringsAsFactors=FALSE)

mohn <- read.csv("results/Mohn_2008_promoters.txt", sep="\t")

histone <- read.csv("results/histone_express_calls.csv")

#convert mohn data to ensembl ID somehow...RefSeq?

#last row doesnt have an id - remove..
histone <- histone[1:15373,]
rownames(histone) <- histone[,"EnsemblID"]

#remove rows that don't have Ref_Seq_IDs - seem to be peaks being called but not annotated??
mohn.peaks <- mohn[which(mohn["Ref_Seq_ID"]!=""),]

ids <- mohn.peaks[,"Ref_Seq_ID"]

library(biomaRt)
ensmart<- useMart("ensembl",dataset="mmusculus_gene_ensembl")
results<- getBM(filters="refseq_dna", values = ids, attributes=c("ensembl_gene_id", "refseq_dna"), mart=ensmart)

#add ensembl ID back to mohn data

mohn.ens <- merge(mohn.peaks, results, by.x = "Ref_Seq_ID", by.y = "refseq_dna")

#remove duplicate hits

mohn.ens <- mohn.ens[order(mohn.ens[,"ensembl_gene_id"], decreasing = TRUE),]
mohn.ens.nodup <- mohn.ens[!duplicated(mohn.ens[,"ensembl_gene_id"]),]

#combine mohn and my data merged by ensembl id

combine <- merge(histone, mohn.ens.nodup, by.x = "EnsemblID", by.y = "ensembl_gene_id")

#tidy up table

keepcols <- c(1,4,5,6,7,8,9,10,11,14,18,35,36,37,38,39,40,41,42,43)

combine <- combine[,keepcols]

write.csv(combine, "results/mohncomparison.csv")

#compare numbers of modifications

NS.K4 <- length(which(!is.na(combine[,"NS_K4"])))

NP_K4 <- length(which(combine[,"K4_call_NP"]=="K4+"))

##find genes that are marked bivalent in each set

NS.bivalent.mine <- intersect(which(!is.na(combine[,"NS_K4"])),which(!is.na(combine[,"NS_K27"])))
N.bivalent.mine <- intersect(which(!is.na(combine[,"N_K4"])),which(!is.na(combine[,"N_K27"])))

NS.bivalent.mine <- combine[NS.bivalent.mine,"EnsemblID"]
N.bivalent.mine <- combine[N.bivalent.mine,"EnsemblID"]

NS.bivalent.mohn <- intersect(which(combine[,"K4_call_NP"]=="K4+"),which(combine[,"K27_call_NP"]=="K27+"))
N.bivalent.mohn <- intersect(which(combine[,"K4_call_TN"]=="K4+"),which(combine[,"K27_call_TN"]=="K27+"))

NS.bivalent.mohn <- combine[NS.bivalent.mohn, "EnsemblID"]
N.bivalent.mohn <- combine[N.bivalent.mohn, "EnsemblID"]

NS.bivalent.both <- intersect(NS.bivalent.mine,NS.bivalent.mohn)
N.bivalent.both <- intersect(N.bivalent.mine,N.bivalent.mohn)

NS.bivalent.dif <- which(!(NS.bivalent.mohn %in% NS.bivalent.mine))
N.bivalent.dif <- which(!(N.bivalent.mohn %in% N.bivalent.mine))

NS.bivalent.dif <- NS.bivalent.mohn[NS.bivalent.dif]
N.bivalent.dif <- N.bivalent.mohn[N.bivalent.dif]

rownames(combine) <- combine[,"EnsemblID"]

NS.bivalent.dif.df <- combine[as.character(NS.bivalent.dif),]

############################################WTF?????????? this only give some of the rows back, rest just list as NA....

######hmm this is quite different....

##compare to mikkelsen data

mik <- read.csv("results/mikkelsen_chromatin_state.csv")

#reannotate to ensembl

mikrefseq <- mik[,"RefSeq.s."]

mik.ens <- getBM(filters="refseq_dna", values = mikrefseq, attributes=c("ensembl_gene_id", "refseq_dna"), mart=ensmart, uniqueRows = TRUE)

##for the time being remove duplicated ensembl ids

mik.ens.or <- mik.ens[order(mik.ens[,"ensembl_gene_id"],decreasing=FALSE),]

mik.ens.or.nodup <- mik.ens.or[!duplicated(mik.ens.or[,"ensembl_gene_id"]),]

#merge ensembl back to mik

mik.res <- merge(mik, mik.ens.or.nodup, by.x = "RefSeq.s.", by.y = "refseq_dna")

##then merge to histone data

combine.mik <- merge(histone, mik.res, by.x = "EnsemblID", by.y = "ensembl_gene_id")

##work out differences

NS.bivalent.mine.mik <- intersect(which(!is.na(combine.mik[,"NS_K4"])),which(!is.na(combine.mik[,"NS_K27"])))

bivalent.mik <- which(combine.mik[,"NPC"]=="K4+K27")

#find genes that are bivalent in both sets

bivalent.mik.dif <- intersect(intersect(which(!is.na(combine.mik[,"NS_K4"])),which(!is.na(combine.mik[,"NS_K27"]))),which(combine.mik[,"NPC"]=="K4+K27"))



##########what are the genes that are differently bivalent in mohn data vs mine?

write.csv(NS.bivalent.dif, "results/NS.bivalent.mohnvsmine.mine_unique.csv")

#get genes that are uniquely bivalent in mohn dataset

NS.bivalent.mohnunique <- NS.bivalent.mohn[which(!(NS.bivalent.mohn %in% NS.bivalent.both))]

NS.bivalent.mohnunique.df <- combine[NS.bivalent.mohnunique,]

################################################################
##################################################################
###try again with ensembl transcript ids then take out ones that arent bivalent

##start from mohn - take out ones that dont map to ensembl transcript, then map them to a gene then remove ones that arent bivalent and blah....

mohn.t <- mohn[which(mohn[,"EnsID"]==""),]







