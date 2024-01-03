rm(list = ls())
source("functions.R")
source("networkFunctions.R")
networks<-GOSimPPI()
load("../data/ChuDataset.RData")

Chu3 <- log2(Chu3+1.1)
integ.l <- DoIntegPPI(exp.m = Chu3, ppiA.m = networks$PPI)

maxSR<-CompMaxSR(integ.l)
result<-CompSRanaPRL(1, exp.m = integ.l$expMC, adj.m = integ.l$adjMC, local = FALSE, maxSR = maxSR)
###Weighted with GO similarity
integ.l$adjMC<-integ.l$adjMC * (1+networks$GOSim[rownames(integ.l$adjMC), colnames(integ.l$adjMC)])/2


head(integ.l$adjMC[,1:5])
sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 1)  #4)
sr.o.normal <- CompSRana(integ.l, local = TRUE, mc.cores = 1)  #4)
save(sr.o, sr.o.normal, file = "Chu3_results_PPIWithGOSimilarity.RData")

