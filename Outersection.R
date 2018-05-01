outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
  # Copy a data.frame to clipboard
  write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
}


PYC <-read.csv("PYC_Age.csv")
dim(PYC)
Age <- (filter(PYC,PYC$pvalue < 0.05))$Gene.name
length(Age)
######
length(zA <- (filter(S,S$PYC_Zanxietyp.value < 0.05))$GeneName) #Significant zAnxiety and its length
length(os <- outersect(Age, zA)) #Age and Anxiety (Union - Intersection)
length(AgeEqualsZa <- intersect(Age,zA)) #Age and Anxiety intesection
length(fromzA <- intersect(os,zA)) # Significant Anxiety without Age effect
zAmiunsAge <- S[is.element(S$GeneName, fromzA),] # filter terms significant in Anxiety but not age
df <- data_frame(zAmiunsAge$GeneName,zAmiunsAge$PYC_Zanxietyestimate,zAmiunsAge$PYC_Zanxietyp.value) #make dataframe of the filtered terms
colnames(df) <- c("GeneName","PYC_Zanxietyestimate","PYC_Zanxietyp.value") # Rename the columns of the dataframe
(df <- df[order(df$PYC_Zanxietyestimate),]) # sort the datframe based on correlation R value
dfN <- (filter(df,df$PYC_Zanxietyestimate < 0))$GeneName #Filter the negatively correlated (Up in Young and Down in the old)
dfP <- (filter(df,df$PYC_Zanxietyestimate > 0))$GeneName #Filter the positively correlated (Down in Young and Up in the old)
write.csv(dfN,file = "zAminusAgeUp.csv") #Save gene list for GO
write.csv(dfP,file = "zAminusAgeDown.csv") #Save gene list for GO
############ List of genes in the intersection###########
dim(ListIntersect <- S[is.element(S$GeneName, AgeEqualsZa),])
ds <- data_frame(ListIntersect$GeneName,ListIntersect$PYC_Zanxietyestimate,ListIntersect$PYC_Zanxietyp.value)
colnames(ds) <- c("GeneName","PYC_Zanxietyestimate","PYC_Zanxietyp.value")
(ds <- ds[order(ds$PYC_Zanxietyestimate),])
dsN <- (filter(ds,ds$PYC_Zanxietyestimate < 0))$GeneName #Filter the negatively correlated (Up in Young and Down in the old)
dsP <- (filter(ds,ds$PYC_Zanxietyestimate > 0))$GeneName #Filter the negatively correlated (Up in Young and Down in the old)
write.csv(dfN,file = "AgeEqualszAUp.csv") #Save gene list for GO
write.csv(dfP,file = "AgeEqualszADown.csv") #Save gene list for GO


#####
zC <- (filter(S,S$PYC_Zcognitionp.value < 0.05))$GeneName
length(os <- outersect(Age, zC))
length(fromzC <- intersect(os,zC))
filtered <- S[is.element(S$GeneName, fromzC),]
df <- data_frame(filtered$GeneName,filtered$PYC_Zcognitionestimate,filtered$PYC_Zcognitionp.value)
colnames(df) <- c("GeneName","PYC_Zcognitionestimate","PYC_Zcognitionp.value")
(df <- df[order(df$PYC_Zcognitionestimate),])
###########MakePlots#####
eset <- load("PYCeset.RData")
eset
bm = pyc
pdata=pData(bm)
edata=exprs(bm)
fdata = fData(bm)

y2 <- removeBatchEffect(edata, pdata$Age)

geneOfInterest <- "Pygm"
#forGGplot <- tbl_df(cbind(edata[geneOfInterest,] , pdata$PC1, pdata$ZAnxity, pdata$ZCognition, pdata$Age))
forGGplot <- tbl_df(cbind(y2[geneOfInterest,] , pdata$PC1, pdata$ZAnxity, pdata$ZCognition, pdata$Age))
colnames(forGGplot) <- c("GeneExp", "PC1", "ZAnxity", "ZCognition","Age")
#forGGplot[forGGplot$Age == "1", "Age"] <- "Young"
#forGGplot[forGGplot$Age == "2", "Age"] <- "Old"

cor.test(x = forGGplot$GeneExp, y = forGGplot$ZAnxity ,alternative = "two.sided",method =  "spearman")

p <- ggplot(data = forGGplot, aes(y=GeneExp, x=ZAnxity)) + 
  geom_point(aes(color=factor(Age))) + geom_smooth(method="lm") + ylab(geneOfInterest) +
  annotate("text", x=Inf, y=-Inf, label=paste("r = ",round(cor(forGGplot$GeneExp, forGGplot$ZAnxity, method = "spearman"), digits=2)), hjust=1, vjust=-1)

p

# quickplot(pdata$ZCognition,edata["Akirin1",], col=1) 
# abline(lm(edata["Akirin1",] ~ pdata$ZCognition))
##############PV#############
PV <-read.csv("PV_Age.csv")
dim(PV)
Age <- (filter(PV,PV$pvalue < 0.05))$Gene.name
######
zA <- (filter(S,S$PV_Zanxietyp.value < 0.05))$GeneName
length(os <- outersect(Age, zA))
length(fromzA <- intersect(os,zA))
filtered <- S[is.element(S$GeneName, fromzA),]
df <- data_frame(filtered$GeneName,filtered$PV_Zanxietyestimate,filtered$PV_Zanxietyp.value)
colnames(df) <- c("GeneName","PV_Zanxietyestimate","PV_Zanxietyp.value")
(df <- df[order(df$PV_Zanxietyestimate),])
#####
zC <- (filter(S,S$PV_Zcognitionp.value < 0.05))$GeneName
length(os <- outersect(Age, zC))
length(fromzC <- intersect(os,zC))
filtered <- S[is.element(S$GeneName, fromzC),]
df <- data_frame(filtered$GeneName,filtered$PV_Zcognitionestimate,filtered$PV_Zcognitionp.value)
colnames(df) <- c("GeneName","PV_Zcognitionestimate","PV_Zcognitionp.value")
(df <- df[order(df$PV_Zcognitionestimate),])
###########SST#########
SST <-read.csv("SST_Age.csv")
dim(SST)
######
zA <- (filter(S,S$SST_Zanxietyp.value < 0.05))$GeneName
length(os <- outersect(Age, zA))
length(fromzA <- intersect(os,zA))
filtered <- S[is.element(S$GeneName, fromzA),]
df <- data_frame(filtered$GeneName,filtered$SST_Zanxietyestimate,filtered$SST_Zanxietyp.value)
colnames(df) <- c("GeneName","SST_Zanxietyestimate","SST_Zanxietyp.value")
(df <- df[order(df$SST_Zanxietyestimate),])
#####
zC <- (filter(S,S$SST_Zcognitionp.value < 0.05))$GeneName
length(os <- outersect(Age, zC))
length(fromzC <- intersect(os,zC))
filtered <- S[is.element(S$GeneName, fromzC),]
df <- data_frame(filtered$GeneName,filtered$SST_Zcognitionestimate,filtered$SST_Zcognitionp.value)
colnames(df) <- c("GeneName","SST_Zcognitionestimate","SST_Zcognitionp.value")
(df <- df[order(df$SST_Zcognitionestimate),])
###########VIP#########
VIP <-read.csv("VIP_Age.csv")
dim(VIP)
Age <- (filter(VIP,VIP$pvalue < 0.05))$Gene.name
######
zA <- (filter(S,S$VIP_Zanxietyp.value < 0.05))$GeneName
length(os <- outersect(Age, zA))
length(fromzA <- intersect(os,zA))
filtered <- S[is.element(S$GeneName, fromzA),]
df <- data_frame(filtered$GeneName,filtered$VIP_Zanxietyestimate,filtered$VIP_Zanxietyp.value)
colnames(df) <- c("GeneName","VIP_Zanxietyestimate","VIP_Zanxietyp.value")
(df <- df[order(df$VIP_Zanxietyestimate),])
#####
zC <- (filter(S,S$VIP_Zcognitionp.value < 0.05))$GeneName
length(os <- outersect(Age, zC))
length(fromzC <- intersect(os,zC))
filtered <- S[is.element(S$GeneName, fromzC),]
df <- data_frame(filtered$GeneName,filtered$VIP_Zcognitionestimate,filtered$VIP_Zcognitionp.value)
colnames(df) <- c("GeneName","VIP_Zcognitionestimate","VIP_Zcognitionp.value")
(df <- df[order(df$VIP_Zcognitionestimate),])

y <- matrix(rnorm(10*9),10,9)
y[,1:3] <- y[,1:3] + 5
batch <- c("A","A","A","B","B","B","C","C","C")
y2 <- removeBatchEffect(edata, pdata$Age)
par(mfrow=c(1,2))
boxplot(as.data.frame(edata),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")

