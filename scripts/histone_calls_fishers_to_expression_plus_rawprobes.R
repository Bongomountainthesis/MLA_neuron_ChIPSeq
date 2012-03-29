#!/usr/local/bin/Rscript

##############take histone mod subsets and see if they are enriched for genes that change expression eg K4 with up, K27 with down and bivalent resolutions

hist.exp <- read.csv("hist_exp_plus_rawexpvalues.csv")

bivalent.NS.genes <- read.csv("bivalent_NS_genes.csv")
bivalent.NS.genes.length <- length(bivalent.NS.genes)

bivalent.N.genes <- read.csv("bivalent_N_genes.csv")
bivalent.N.genes.length <- length(bivalent.N.genes)

bivalent.NS.kept.genes <- read.csv("bivalent_NS_kept_genes.csv")

bivalent.NS.lost.genes <- read.csv("bivalent_NS_lost_genes.csv")

bivalent.N.gain.genes <- read.csv("bivalent_N_gain_genes.csv")

bivalent.NS.lost.K4.genes <- read.csv("bivalent_NS_lost_to_K4_genes.csv")

bivalent.NS.lost.K27.genes <- read.csv("bivalent_NS_lost_K27_genes.csv")

bivalent.NS.lost.nothing.genes <- read.csv("bivalent_NS_lost_nothing_genes.csv")

################also calculate genes that just have K4 or K27 in NS and N

NS.K4 <- intersect(which(hist.exp[,"NS_K4"]!="<NA>"), which(is.na(hist.exp[,"NS_K27"])))
NS.K27 <- intersect(which(hist.exp[,"NS_K27"]!="<NA>"), which(is.na(hist.exp[,"NS_K4"])))

N.K4 <- intersect(which(hist.exp[,"N_K4"]!="<NA>"), which(is.na(hist.exp[,"N_K27"])))
N.K27 <- intersect(which(hist.exp[,"N_K27"]!="<NA>"), which(is.na(hist.exp[,"N_K4"])))

NS.K4 <- hist.exp[NS.K4,]
NS.K27 <- hist.exp[NS.K27,]

N.K4 <- hist.exp[N.K4,]
N.K27 <- hist.exp[N.K27,]

NS.K4.length <- length(NS.K4[,"EnsemblID"])
NS.K27.length <- length(NS.K27[,"EnsemblID"])

N.K4.length <- length(N.K4[,"EnsemblID"])
N.K27.length <- length(N.K27[,"EnsemblID"])

##################subset genes into up and down regulated

##all genes into up and down

up.genes <- hist.exp[which(hist.exp[,"logFC"]>=1),]
down.genes <- hist.exp[which(hist.exp[,"logFC"]<=-1),]

up <- length(up.genes[,"EnsemblID"])
down <- length(down.genes[,"EnsemblID"])
all <- length(hist.exp[,"EnsemblID"])

##all K4 and K27 up and down

NS.K4.up <- intersect(NS.K4[,"EnsemblID"], up.genes[,"EnsemblID"])
NS.K4.down <- intersect(NS.K4[,"EnsemblID"], down.genes[,"EnsemblID"])

NS.K27.up <- intersect(NS.K27[,"EnsemblID"], up.genes[,"EnsemblID"])
NS.K27.down <- intersect(NS.K27[,"EnsemblID"], down.genes[,"EnsemblID"])

NS.K4.up.length <- length(NS.K4.up)
NS.K4.down.length <- length(NS.K4.down)

NS.K27.up.length <- length(NS.K27.up)
NS.K27.down.length <- length(NS.K27.down)

N.K4.up <- intersect(N.K4[,"EnsemblID"], up.genes[,"EnsemblID"])
N.K4.down <- intersect(N.K4[,"EnsemblID"], down.genes[,"EnsemblID"])

N.K27.up <- intersect(N.K27[,"EnsemblID"], up.genes[,"EnsemblID"])
N.K27.down <- intersect(N.K27[,"EnsemblID"], down.genes[,"EnsemblID"])

N.K4.up.length <- length(N.K4.up)
N.K4.down.length <- length(N.K4.down)

N.K27.up.length <- length(N.K27.up)
N.K27.down.length <- length(N.K27.down)


##############then split subsets by genes that are up/down regulated
##get genes

bivalent.NS.kept.upl <- intersect(bivalent.NS.kept.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.NS.kept.downl <- intersect(bivalent.NS.kept.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.N.gain.upl <- intersect(bivalent.N.gain.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.N.gain.downl <- intersect(bivalent.N.gain.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.NS.lost.upl <- intersect(bivalent.NS.lost.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.NS.lost.downl <- intersect(bivalent.NS.lost.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.NS.lost.K4.upl <- intersect(bivalent.NS.lost.K4.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.NS.lost.K4.downl <- intersect(bivalent.NS.lost.K4.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.NS.lost.K27.upl <- intersect(bivalent.NS.lost.K27.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.NS.lost.K27.downl <- intersect(bivalent.NS.lost.K27.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.NS.lost.nothing.upl <- intersect(bivalent.NS.lost.nothing.genes[,"EnsemblID"], up.genes[,"EnsemblID"])
bivalent.NS.lost.nothing.downl <- intersect(bivalent.NS.lost.nothing.genes[,"EnsemblID"], down.genes[,"EnsemblID"])

bivalent.NS.kept.up <- length(bivalent.NS.kept.upl)
bivalent.NS.kept.down <- length(bivalent.NS.kept.downl)

bivalent.N.gain.up <- length(bivalent.N.gain.upl)
bivalent.N.gain.down <- length(bivalent.N.gain.downl)

bivalent.NS.lost.up <- length(bivalent.NS.lost.upl)
bivalent.NS.lost.down <- length(bivalent.NS.lost.downl)

bivalent.NS.lost.K4.up <- length(bivalent.NS.lost.K4.upl)
bivalent.NS.lost.K4.down <- length(bivalent.NS.lost.K4.downl)

bivalent.NS.lost.K27.up <- length(bivalent.NS.lost.K27.upl)
bivalent.NS.lost.K27.down <- length(bivalent.NS.lost.K27.downl)

bivalent.NS.lost.nothing.up <- length(bivalent.NS.lost.nothing.upl)
bivalent.NS.lost.nothing.down <- length(bivalent.NS.lost.nothing.downl)

###find genes are/arent expressed

NS.notexpressed <- hist.exp[which(hist.exp[,"NS"]<=8),"EnsemblID"]

N.notexpressed <- hist.exp[which(hist.exp[,"N"]<=8),"EnsemblID"]

NS.expressed <- hist.exp[which(hist.exp[,"NS"]>=8),"EnsemblID"]

N.expressed <- hist.exp[which(hist.exp[,"N"]>=8), "EnsemblID"]

####link to histones
#genes that have K4 and are expressed

NS.K4.expressed <- which(NS.expressed %in% NS.K4[,"EnsemblID"])
N.K4.expressed <- which(N.expressed %in% N.K4[,"EnsemblID"])

NS.noK4.expressed <- NS.expressed[which(!(NS.expressed %in% NS.K4[,"EnsemblID"]))]
NS.noK4.expressed <- hist.exp[hist.exp[,"EnsemblID"] %in% NS.noK4.expressed,]

NS.noK4.expressed.bivalent <- which(NS.noK4.expressed %in% bivalent.NS.genes[,"EnsemblID"])
NS.noK4.expressed.K27 <- which(NS.noK4.expressed %in% NS.K27[,"EnsemblID"])

N.noK4.expressed <- N.expressed[which(!(N.expressed %in% N.K4[,"EnsemblID"]))]
N.noK4.expressed <- hist.exp[hist.exp[,"EnsemblID"] %in% N.noK4.expressed,]

N.noK4.expressed.bivalent <- which(N.noK4.expressed %in% bivalent.N.genes[,"EnsemblID"])
N.noK4.expressed.K27 <- which(N.noK4.expressed %in% N.K27[,"EnsemblID"])

NS.K27.notexpressed <- which(NS.notexpressed %in% NS.K27[,"EnsemblID"])
N.K27.notexpressed <- which(N.notexpressed %in% N.K27[,"EnsemblID"])

NS.K27.expressed <- which(NS.expressed %in% NS.K27[,"EnsemblID"])
N.K27.expressed <- which(N.expressed %in% N.K27[,"EnsemblID"])

NS.bivalent.expressed <- which(NS.expressed %in% bivalent.NS.genes[,"EnsemblID"])
N.bivalent.expressed <- which(N.expressed %in% bivalent.N.genes[,"EnsemblID"])

NS.bivalent.notexpressed <- which(NS.notexpressed %in% bivalent.NS.genes[,"EnsemblID"])
N.bivalent.notexpressed <- which(N.notexpressed %in% bivalent.N.genes[,"EnsemblID"])



################make contingency matrix

hist.up <- bivalent.NS.lost.K4.down
hist.total <- bivalent.NS.genes.length

pop.up <- up
pop.total <- all

matrix <- matrix(c(hist.up, hist.total, pop.up, pop.total), nr = 2)

matrix <- matrix(c(, 634, down, all), nr = 2)

fisher.test(matrix)














