pvA <- read.csv("PVAxnietyAgeC.csv")
pvC <- read.csv("PVCognitionAgeC.csv")
removeFromCorrelationPV <- intersect(pvC$GeneName,pvA$GeneName)
SremovedPVIntersection <- S[-which(S$GeneName %in% removeFromCorrelationPV), ]
cor.test(S$PV_Zanxietyestimate,S$PV_Zcognitionestimate,method = "spearman")
cor.test(SremovedPVIntersection$PV_Zanxietyestimate,SremovedPVIntersection$PV_Zcognitionestimate, method = "spearman")

#read Age data 
PYC <-read.csv("PYC_Age.csv")
dim(PYC)
L <- (filter(PYC,PYC$pvalue < 0.05))$Gene.name
U <- (filter(PYC,PYC$pvalue < 0.05,PYC$stat > 0))$Gene.name
D <- (filter(PYC,PYC$pvalue < 0.05,PYC$stat < 0))$Gene.name
length(U) + length(D)
######
pyc1 <- (filter(S,S$PYC_PC1p.value < 0.05 ))$GeneName
pyc1.u0.05 <- (filter(S,S$PYC_PC1p.value < 0.05 & S$PYC_PC1estimate > 0))$GeneName
pyc1.u0.01 <- (filter(S,S$PYC_PC1p.value < 0.01 & S$PYC_PC1estimate > 0))$GeneName
length(pyc1.u0.05);length(pyc1.u0.01)
pyc1.d0.05 <- (filter(S,S$PYC_PC1p.value < 0.05 & S$PYC_PC1estimate < 0))$GeneName
pyc1.d0.01 <- (filter(S,S$PYC_PC1p.value < 0.01 & S$PYC_PC1estimate < 0))$GeneName
length(pyc1.d0.05); length(pyc1.d0.01);  length(pyc1.u0.05); length(pyc1.u0.01)
length(pyc1.d0.05) + length(pyc1.u0.05); length(pyc1.d0.01) + length(pyc1.u0.01)
length(intersect(pyc1.u0.05,D));length(intersect(pyc1.u0.01,D));length(intersect(pyc1.d0.05,U));length(intersect(pyc1.d0.01,U))
length(intersect(pyc1.u0.05,U));length(intersect(pyc1.u0.01,U));length(intersect(pyc1.d0.05,D));length(intersect(pyc1.d0.01,D))
dim(pyc1)
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

length(os <- outersect(L, pyc1))
length(frompyc <- intersect(os,pyc1))
filtered <- S[is.element(S$GeneName, frompyc),]
df <- data_frame(filtered$GeneName,filtered$PYC_PC1estimate,filtered$PYC_PC1p.value)
df <- df[order(df$`filtered$PYC_PC1estimate`),] 
#(S[,S$PYC_PC1p.value %in% frompyc])
#filter(df, S %in% frompyc)









