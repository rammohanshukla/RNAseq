P <- read.csv("Pearson.csv")
S <- read.csv("Spearman.csv")
cor.test(P$PV_PC1estimate,P$PV_PC2estimate)
cor.test(P$PV_PC1estimate,P$PV_PC3estimate)
cor.test(P$PV_PC1estimate,P$PV_Zanxietyestimate)
cor.test(P$PV_PC1estimate,P$PV_Zcognitionestimate)
cor.test(P$PV_PC2estimate,P$PV_PC3estimate)
cor.test(P$PV_PC2estimate,P$PV_Zanxietyestimate)
cor.test(P$PV_PC2estimate,P$PV_Zcognitionestimate)
cor.test(P$PV_PC3estimate,P$PV_Zanxietyestimate)
cor.test(P$PV_PC3estimate,P$PV_Zcognitionestimate)

pvA <- read.csv("PVAxnietyAgeC.csv")
pvC <- read.csv("PVCognitionAgeC.csv")
removeFromCorrelationPV <- intersect(pvC$GeneName,pvA$GeneName)
SremovedPVIntersection <- S[-which(S$GeneName %in% removeFromCorrelationPV), ]
cor.test(S$PV_Zanxietyestimate,S$PV_Zcognitionestimate,method = "spearman")
cor.test(SremovedPVIntersection$PV_Zanxietyestimate,SremovedPVIntersection$PV_Zcognitionestimate, method = "spearman")


cor.test(P$PV_PC1estimate,S$PV_PC1estimate)
cor.test(P$PV_PC2estimate,S$PV_PC2estimate)
cor.test(P$PV_PC3estimate,S$PV_PC3estimate)
cor.test(P$PV_Zanxietyestimate,S$PV_Zanxietyestimate)
cor.test(P$PV_Zcognitionestimate,S$PV_Zcognitionestimate)



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

######
pyc2.u0.05 <- (filter(S,S$PYC_PC2p.value < 0.05 & S$PYC_PC2estimate > 0))$GeneName
pyc2.u0.01 <- (filter(S,S$PYC_PC2p.value < 0.01 & S$PYC_PC2estimate > 0))$GeneName
length(pyc2.u0.05);length(pyc2.u0.01)
pyc2.d0.05 <- (filter(S,S$PYC_PC2p.value < 0.05 & S$PYC_PC2estimate < 0))$GeneName
pyc2.d0.01 <- (filter(S,S$PYC_PC2p.value < 0.01 & S$PYC_PC2estimate < 0))$GeneName
length(pyc2.d0.05); length(pyc2.d0.01);  length(pyc2.u0.05); length(pyc2.u0.01)
length(pyc2.d0.05) + length(pyc2.u0.05); length(pyc2.d0.01) + length(pyc2.u0.01)
length(intersect(pyc2.u0.05,D));length(intersect(pyc2.u0.01,D));length(intersect(pyc2.d0.05,U));length(intersect(pyc2.d0.01,U))
length(intersect(pyc2.u0.05,U));length(intersect(pyc2.u0.01,U));length(intersect(pyc2.d0.05,D));length(intersect(pyc2.d0.01,D))
######
pyc3.u0.05 <- (filter(S,S$PYC_PC3p.value < 0.05 & S$PYC_PC3estimate > 0))$GeneName
pyc3.u0.01 <- (filter(S,S$PYC_PC3p.value < 0.01 & S$PYC_PC3estimate > 0))$GeneName
length(pyc3.u0.05);length(pyc3.u0.01)
pyc3.d0.05 <- (filter(S,S$PYC_PC3p.value < 0.05 & S$PYC_PC3estimate < 0))$GeneName
pyc3.d0.01 <- (filter(S,S$PYC_PC3p.value < 0.01 & S$PYC_PC3estimate < 0))$GeneName
length(pyc3.d0.05); length(pyc3.d0.01);  length(pyc3.u0.05); length(pyc3.u0.01)
length(pyc3.d0.05) + length(pyc3.u0.05); length(pyc3.d0.01) + length(pyc3.u0.01)
length(intersect(pyc3.u0.05,D));length(intersect(pyc3.u0.01,D));length(intersect(pyc3.d0.05,U));length(intersect(pyc3.d0.01,U))
length(intersect(pyc3.u0.05,U));length(intersect(pyc3.u0.01,U));length(intersect(pyc3.d0.05,D));length(intersect(pyc3.d0.01,D))
######
zA.u0.05 <- (filter(S,S$PYC_Zanxietyp.value < 0.05 & S$PYC_Zanxietyestimate > 0))$GeneName
zA.u0.01 <- (filter(S,S$PYC_Zanxietyp.value < 0.01 & S$PYC_Zanxietyestimate > 0))$GeneName
length(zA.u0.05);length(zA.u0.01)
zA.d0.05 <- (filter(S,S$PYC_Zanxietyp.value < 0.05 & S$PYC_Zanxietyestimate < 0))$GeneName
zA.d0.01 <- (filter(S,S$PYC_Zanxietyp.value < 0.01 & S$PYC_Zanxietyestimate < 0))$GeneName
length(zA.d0.05); length(zA.d0.01);  length(zA.u0.05); length(zA.u0.01)
length(zA.d0.05) + length(zA.u0.05); length(zA.d0.01) + length(zA.u0.01)
length(intersect(zA.u0.05,D));length(intersect(zA.u0.01,D));length(intersect(zA.d0.05,U));length(intersect(zA.d0.01,U))
length(intersect(zA.u0.05,U));length(intersect(zA.u0.01,U));length(intersect(zA.d0.05,D));length(intersect(zA.d0.01,D))
#####
zC.u0.05 <- (filter(S,S$PYC_Zcognitionp.value < 0.05 & S$PYC_Zcognitionestimate > 0))$GeneName
zC.u0.01 <- (filter(S,S$PYC_Zcognitionp.value < 0.01 & S$PYC_Zcognitionestimate > 0))$GeneName
length(zC.u0.05);length(zC.u0.01)
zC.d0.05 <- (filter(S,S$PYC_Zcognitionp.value < 0.05 & S$PYC_Zcognitionestimate < 0))$GeneName
zC.d0.01 <- (filter(S,S$PYC_Zcognitionp.value < 0.01 & S$PYC_Zcognitionestimate < 0))$GeneName
length(zC.d0.05); length(zC.d0.01);  length(zC.u0.05); length(zC.u0.01)
length(zC.d0.05) + length(zC.u0.05); length(zC.d0.01) + length(zC.u0.01)
length(intersect(zC.u0.05,D));length(intersect(zC.u0.01,D));length(intersect(zC.d0.05,U));length(intersect(zC.d0.01,U))
length(intersect(zC.u0.05,U));length(intersect(zC.u0.01,U));length(intersect(zC.d0.05,D));length(intersect(zC.d0.01,D))
##########################PV########################
PV <-read.csv("PV_Age.csv")
dim(PV)
U <- (filter(PV,PV$pvalue < 0.05,PV$stat > 0))$Gene.name
D <- (filter(PV,PV$pvalue < 0.05,PV$stat < 0))$Gene.name
length(U1) + length(D1)
######
pv1.u0.05 <- (filter(S,S$PV_PC1p.value < 0.05 & S$PV_PC1estimate > 0))$GeneName
pv1.u0.01 <- (filter(S,S$PV_PC1p.value < 0.01 & S$PV_PC1estimate > 0))$GeneName
length(pv1.u0.05);length(pv1.u0.01)
pv1.d0.05 <- (filter(S,S$PV_PC1p.value < 0.05 & S$PV_PC1estimate < 0))$GeneName
pv1.d0.01 <- (filter(S,S$PV_PC1p.value < 0.01 & S$PV_PC1estimate < 0))$GeneName
length(pv1.d0.05); length(pv1.d0.01);  length(pv1.u0.05); length(pv1.u0.01)
length(pv1.d0.05) + length(pv1.u0.05); length(pv1.d0.01) + length(pv1.u0.01)
length(intersect(pv1.u0.05,D));length(intersect(pv1.u0.01,D));length(intersect(pv1.d0.05,U));length(intersect(pv1.d0.01,U))
length(intersect(pv1.u0.05,U));length(intersect(pv1.u0.01,U));length(intersect(pv1.d0.05,D));length(intersect(pv1.d0.01,D))
######
pv2.u0.05 <- (filter(S,S$PV_PC2p.value < 0.05 & S$PV_PC2estimate > 0))$GeneName
pv2.u0.01 <- (filter(S,S$PV_PC2p.value < 0.01 & S$PV_PC2estimate > 0))$GeneName
length(pv2.u0.05);length(pv2.u0.01)
pv2.d0.05 <- (filter(S,S$PV_PC2p.value < 0.05 & S$PV_PC2estimate < 0))$GeneName
pv2.d0.01 <- (filter(S,S$PV_PC2p.value < 0.01 & S$PV_PC2estimate < 0))$GeneName
length(pv2.d0.05); length(pv2.d0.01);  length(pv2.u0.05); length(pv2.u0.01)
length(pv2.d0.05) + length(pv2.u0.05); length(pv2.d0.01) + length(pv2.u0.01)
length(intersect(pv2.u0.05,D));length(intersect(pv2.u0.01,D));length(intersect(pv2.d0.05,U));length(intersect(pv2.d0.01,U))
length(intersect(pv2.u0.05,U));length(intersect(pv2.u0.01,U));length(intersect(pv2.d0.05,D));length(intersect(pv2.d0.01,D))
######
pv3.u0.05 <- (filter(S,S$PV_PC3p.value < 0.05 & S$PV_PC3estimate > 0))$GeneName
pv3.u0.01 <- (filter(S,S$PV_PC3p.value < 0.01 & S$PV_PC3estimate > 0))$GeneName
length(pv3.u0.05);length(pv3.u0.01)
pv3.d0.05 <- (filter(S,S$PV_PC3p.value < 0.05 & S$PV_PC3estimate < 0))$GeneName
pv3.d0.01 <- (filter(S,S$PV_PC3p.value < 0.01 & S$PV_PC3estimate < 0))$GeneName
length(pv3.d0.05); length(pv3.d0.01);  length(pv3.u0.05); length(pv3.u0.01)
length(pv3.d0.05) + length(pv3.u0.05); length(pv3.d0.01) + length(pv3.u0.01)
length(intersect(pv3.u0.05,D));length(intersect(pv3.u0.01,D));length(intersect(pv3.d0.05,U));length(intersect(pv3.d0.01,U))
length(intersect(pv3.u0.05,U));length(intersect(pv3.u0.01,U));length(intersect(pv3.d0.05,D));length(intersect(pv3.d0.01,D))
######
zA.u0.05 <- (filter(S,S$PV_Zanxietyp.value < 0.05 & S$PV_Zanxietyestimate > 0))$GeneName
zA.u0.01 <- (filter(S,S$PV_Zanxietyp.value < 0.01 & S$PV_Zanxietyestimate > 0))$GeneName
length(zA.u0.05);length(zA.u0.01)
zA.d0.05 <- (filter(S,S$PV_Zanxietyp.value < 0.05 & S$PV_Zanxietyestimate < 0))$GeneName
zA.d0.01 <- (filter(S,S$PV_Zanxietyp.value < 0.01 & S$PV_Zanxietyestimate < 0))$GeneName
length(zA.d0.05); length(zA.d0.01);  length(zA.u0.05); length(zA.u0.01)
length(zA.d0.05) + length(zA.u0.05); length(zA.d0.01) + length(zA.u0.01)
length(intersect(zA.u0.05,D));length(intersect(zA.u0.01,D));length(intersect(zA.d0.05,U));length(intersect(zA.d0.01,U))
length(intersect(zA.u0.05,U));length(intersect(zA.u0.01,U));length(intersect(zA.d0.05,D));length(intersect(zA.d0.01,D))
#####
zC.u0.05 <- (filter(S,S$PV_Zcognitionp.value < 0.05 & S$PV_Zcognitionestimate > 0))$GeneName
zC.u0.01 <- (filter(S,S$PV_Zcognitionp.value < 0.01 & S$PV_Zcognitionestimate > 0))$GeneName
length(zC.u0.05);length(zC.u0.01)
zC.d0.05 <- (filter(S,S$PV_Zcognitionp.value < 0.05 & S$PV_Zcognitionestimate < 0))$GeneName
zC.d0.01 <- (filter(S,S$PV_Zcognitionp.value < 0.01 & S$PV_Zcognitionestimate < 0))$GeneName
length(zC.d0.05); length(zC.d0.01);  length(zC.u0.05); length(zC.u0.01)
length(zC.d0.05) + length(zC.u0.05); length(zC.d0.01) + length(zC.u0.01)
length(intersect(zC.u0.05,D));length(intersect(zC.u0.01,D));length(intersect(zC.d0.05,U));length(intersect(zC.d0.01,U))
length(intersect(zC.u0.05,U));length(intersect(zC.u0.01,U));length(intersect(zC.d0.05,D));length(intersect(zC.d0.01,D))
########################SST#########################
SST <-read.csv("SST_Age.csv")
dim(SST)
U <- (filter(SST,SST$pvalue < 0.05,SST$stat > 0))$Gene.name
D <- (filter(SST,SST$pvalue < 0.05,SST$stat < 0))$Gene.name
length(U) + length(D)
######
sst1.u0.05 <- (filter(S,S$SST_PC1p.value < 0.05 & S$SST_PC1estimate > 0))$GeneName
sst1.u0.01 <- (filter(S,S$SST_PC1p.value < 0.01 & S$SST_PC1estimate > 0))$GeneName
length(sst1.u0.05);length(sst1.u0.01)
sst1.d0.05 <- (filter(S,S$SST_PC1p.value < 0.05 & S$SST_PC1estimate < 0))$GeneName
sst1.d0.01 <- (filter(S,S$SST_PC1p.value < 0.01 & S$SST_PC1estimate < 0))$GeneName
length(sst1.d0.05); length(sst1.d0.01);  length(sst1.u0.05); length(sst1.u0.01)
length(sst1.d0.05) + length(sst1.u0.05); length(sst1.d0.01) + length(sst1.u0.01)
length(intersect(sst1.u0.05,D));length(intersect(sst1.u0.01,D));length(intersect(sst1.d0.05,U));length(intersect(sst1.d0.01,U))
length(intersect(sst1.u0.05,U));length(intersect(sst1.u0.01,U));length(intersect(sst1.d0.05,D));length(intersect(sst1.d0.01,D))
######
sst2.u0.05 <- (filter(S,S$SST_PC2p.value < 0.05 & S$SST_PC2estimate > 0))$GeneName
sst2.u0.01 <- (filter(S,S$SST_PC2p.value < 0.01 & S$SST_PC2estimate > 0))$GeneName
length(sst2.u0.05);length(sst2.u0.01)
sst2.d0.05 <- (filter(S,S$SST_PC2p.value < 0.05 & S$SST_PC2estimate < 0))$GeneName
sst2.d0.01 <- (filter(S,S$SST_PC2p.value < 0.01 & S$SST_PC2estimate < 0))$GeneName
length(sst2.d0.05); length(sst2.d0.01);  length(sst2.u0.05); length(sst2.u0.01)
length(sst2.d0.05) + length(sst2.u0.05); length(sst2.d0.01) + length(sst2.u0.01)
length(intersect(sst2.u0.05,D));length(intersect(sst2.u0.01,D));length(intersect(sst2.d0.05,U));length(intersect(sst2.d0.01,U))
length(intersect(sst2.u0.05,U));length(intersect(sst2.u0.01,U));length(intersect(sst2.d0.05,D));length(intersect(sst2.d0.01,D))
######
sst3.u0.05 <- (filter(S,S$SST_PC3p.value < 0.05 & S$SST_PC3estimate > 0))$GeneName
sst3.u0.01 <- (filter(S,S$SST_PC3p.value < 0.01 & S$SST_PC3estimate > 0))$GeneName
length(sst3.u0.05);length(sst3.u0.01)
sst3.d0.05 <- (filter(S,S$SST_PC3p.value < 0.05 & S$SST_PC3estimate < 0))$GeneName
sst3.d0.01 <- (filter(S,S$SST_PC3p.value < 0.01 & S$SST_PC3estimate < 0))$GeneName
length(sst3.d0.05); length(sst3.d0.01);  length(sst3.u0.05); length(sst3.u0.01)
length(sst3.d0.05) + length(sst3.u0.05); length(sst3.d0.01) + length(sst3.u0.01)
length(intersect(sst3.u0.05,D));length(intersect(sst3.u0.01,D));length(intersect(sst3.d0.05,U));length(intersect(sst3.d0.01,U))
length(intersect(sst3.u0.05,U));length(intersect(sst3.u0.01,U));length(intersect(sst3.d0.05,D));length(intersect(sst3.d0.01,D))
######
zA.u0.05 <- (filter(S,S$SST_Zanxietyp.value < 0.05 & S$SST_Zanxietyestimate > 0))$GeneName
zA.u0.01 <- (filter(S,S$SST_Zanxietyp.value < 0.01 & S$SST_Zanxietyestimate > 0))$GeneName
length(zA.u0.05);length(zA.u0.01)
zA.d0.05 <- (filter(S,S$SST_Zanxietyp.value < 0.05 & S$SST_Zanxietyestimate < 0))$GeneName
zA.d0.01 <- (filter(S,S$SST_Zanxietyp.value < 0.01 & S$SST_Zanxietyestimate < 0))$GeneName
length(zA.d0.05); length(zA.d0.01);  length(zA.u0.05); length(zA.u0.01)
length(zA.d0.05) + length(zA.u0.05); length(zA.d0.01) + length(zA.u0.01)
length(intersect(zA.u0.05,D));length(intersect(zA.u0.01,D));length(intersect(zA.d0.05,U));length(intersect(zA.d0.01,U))
length(intersect(zA.u0.05,U));length(intersect(zA.u0.01,U));length(intersect(zA.d0.05,D));length(intersect(zA.d0.01,D))
#####
zC.u0.05 <- (filter(S,S$SST_Zcognitionp.value < 0.05 & S$SST_Zcognitionestimate > 0))$GeneName
zC.u0.01 <- (filter(S,S$SST_Zcognitionp.value < 0.01 & S$SST_Zcognitionestimate > 0))$GeneName
length(zC.u0.05);length(zC.u0.01)
zC.d0.05 <- (filter(S,S$SST_Zcognitionp.value < 0.05 & S$SST_Zcognitionestimate < 0))$GeneName
zC.d0.01 <- (filter(S,S$SST_Zcognitionp.value < 0.01 & S$SST_Zcognitionestimate < 0))$GeneName
length(zC.d0.05); length(zC.d0.01);  length(zC.u0.05); length(zC.u0.01)
length(zC.d0.05) + length(zC.u0.05); length(zC.d0.01) + length(zC.u0.01)
length(intersect(zC.u0.05,D));length(intersect(zC.u0.01,D));length(intersect(zC.d0.05,U));length(intersect(zC.d0.01,U))
length(intersect(zC.u0.05,U));length(intersect(zC.u0.01,U));length(intersect(zC.d0.05,D));length(intersect(zC.d0.01,D))
########################VIP####################
VIP <-read.csv("VIP_Age.csv")
dim(VIP)
U <- (filter(VIP,VIP$pvalue < 0.05,VIP$stat > 0))$Gene.name
D <- (filter(VIP,VIP$pvalue < 0.05,VIP$stat < 0))$Gene.name
length(U) + length(D)
######
vip1.u0.05 <- (filter(S,S$VIP_PC1p.value < 0.05 & S$VIP_PC1estimate > 0))$GeneName
vip1.u0.01 <- (filter(S,S$VIP_PC1p.value < 0.01 & S$VIP_PC1estimate > 0))$GeneName
length(vip1.u0.05);length(vip1.u0.01)
vip1.d0.05 <- (filter(S,S$VIP_PC1p.value < 0.05 & S$VIP_PC1estimate < 0))$GeneName
vip1.d0.01 <- (filter(S,S$VIP_PC1p.value < 0.01 & S$VIP_PC1estimate < 0))$GeneName
length(vip1.d0.05); length(vip1.d0.01);  length(vip1.u0.05); length(vip1.u0.01)
length(vip1.d0.05) + length(vip1.u0.05); length(vip1.d0.01) + length(vip1.u0.01)
length(intersect(vip1.u0.05,D));length(intersect(vip1.u0.01,D));length(intersect(vip1.d0.05,U));length(intersect(vip1.d0.01,U))
length(intersect(vip1.u0.05,U));length(intersect(vip1.u0.01,U));length(intersect(vip1.d0.05,D));length(intersect(vip1.d0.01,D))
######
vip2.u0.05 <- (filter(S,S$VIP_PC2p.value < 0.05 & S$VIP_PC2estimate > 0))$GeneName
vip2.u0.01 <- (filter(S,S$VIP_PC2p.value < 0.01 & S$VIP_PC2estimate > 0))$GeneName
length(vip2.u0.05);length(vip2.u0.01)
vip2.d0.05 <- (filter(S,S$VIP_PC2p.value < 0.05 & S$VIP_PC2estimate < 0))$GeneName
vip2.d0.01 <- (filter(S,S$VIP_PC2p.value < 0.01 & S$VIP_PC2estimate < 0))$GeneName
length(vip2.d0.05); length(vip2.d0.01);  length(vip2.u0.05); length(vip2.u0.01)
length(vip2.d0.05) + length(vip2.u0.05); length(vip2.d0.01) + length(vip2.u0.01)
length(intersect(vip2.u0.05,D));length(intersect(vip2.u0.01,D));length(intersect(vip2.d0.05,U));length(intersect(vip2.d0.01,U))
length(intersect(vip2.u0.05,U));length(intersect(vip2.u0.01,U));length(intersect(vip2.d0.05,D));length(intersect(vip2.d0.01,D))
######
vip3.u0.05 <- (filter(S,S$VIP_PC3p.value < 0.05 & S$VIP_PC3estimate > 0))$GeneName
vip3.u0.01 <- (filter(S,S$VIP_PC3p.value < 0.01 & S$VIP_PC3estimate > 0))$GeneName
length(vip3.u0.05);length(vip3.u0.01)
vip3.d0.05 <- (filter(S,S$VIP_PC3p.value < 0.05 & S$VIP_PC3estimate < 0))$GeneName
vip3.d0.01 <- (filter(S,S$VIP_PC3p.value < 0.01 & S$VIP_PC3estimate < 0))$GeneName
length(vip3.d0.05); length(vip3.d0.01);  length(vip3.u0.05); length(vip3.u0.01)
length(vip3.d0.05) + length(vip3.u0.05); length(vip3.d0.01) + length(vip3.u0.01)
length(intersect(vip3.u0.05,D));length(intersect(vip3.u0.01,D));length(intersect(vip3.d0.05,U));length(intersect(vip3.d0.01,U))
length(intersect(vip3.u0.05,U));length(intersect(vip3.u0.01,U));length(intersect(vip3.d0.05,D));length(intersect(vip3.d0.01,D))
######
zA.u0.05 <- (filter(S,S$VIP_Zanxietyp.value < 0.05 & S$VIP_Zanxietyestimate > 0))$GeneName
zA.u0.01 <- (filter(S,S$VIP_Zanxietyp.value < 0.01 & S$VIP_Zanxietyestimate > 0))$GeneName
length(zA.u0.05);length(zA.u0.01)
zA.d0.05 <- (filter(S,S$VIP_Zanxietyp.value < 0.05 & S$VIP_Zanxietyestimate < 0))$GeneName
zA.d0.01 <- (filter(S,S$VIP_Zanxietyp.value < 0.01 & S$VIP_Zanxietyestimate < 0))$GeneName
length(zA.d0.05); length(zA.d0.01);  length(zA.u0.05); length(zA.u0.01)
length(zA.d0.05) + length(zA.u0.05); length(zA.d0.01) + length(zA.u0.01)
length(intersect(zA.u0.05,D));length(intersect(zA.u0.01,D));length(intersect(zA.d0.05,U));length(intersect(zA.d0.01,U))
length(intersect(zA.u0.05,U));length(intersect(zA.u0.01,U));length(intersect(zA.d0.05,D));length(intersect(zA.d0.01,D))
#####
zC.u0.05 <- (filter(S,S$VIP_Zcognitionp.value < 0.05 & S$VIP_Zcognitionestimate > 0))$GeneName
zC.u0.01 <- (filter(S,S$VIP_Zcognitionp.value < 0.01 & S$VIP_Zcognitionestimate > 0))$GeneName
length(zC.u0.05);length(zC.u0.01)
zC.d0.05 <- (filter(S,S$VIP_Zcognitionp.value < 0.05 & S$VIP_Zcognitionestimate < 0))$GeneName
zC.d0.01 <- (filter(S,S$VIP_Zcognitionp.value < 0.01 & S$VIP_Zcognitionestimate < 0))$GeneName
length(zC.d0.05); length(zC.d0.01);  length(zC.u0.05); length(zC.u0.01)
length(zC.d0.05) + length(zC.u0.05); length(zC.d0.01) + length(zC.u0.01)
length(intersect(zC.u0.05,D));length(intersect(zC.u0.01,D));length(intersect(zC.d0.05,U));length(intersect(zC.d0.01,U))
length(intersect(zC.u0.05,U));length(intersect(zC.u0.01,U));length(intersect(zC.d0.05,D));length(intersect(zC.d0.01,D))


df <- read.csv("epDataPYCGRCm38.csv",header=T)
dft <- (t(df[,-1]))
head(dft[,1:6])
res <- apply(dft[,-c(1:6)], 2, function(x,y) cor.test(x,y,method = "spearman",continuity = FALSE,conf.level = 0.95),dft[,2])
PYC_ZAxniety <- data.frame(do.call(rbind,res))
PYC_ZAxniety <- apply(PYC_ZAxniety,2,as.character)
write.csv(PYC_ZAxniety,file="PYC_ZAxniety.csv")

cor(df[2,],df[4,])




