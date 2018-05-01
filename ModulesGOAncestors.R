library(readr)
library(dplyr)
library(magrittr)
library(GO.db)
library(reshape2)
library(plyr)

list.files(pattern = "*.txt")
#resultGroupsUp <- read_tsv("pathways2experiments_significant_Up_0.05.txt") %>% rename(X1 = "name")
#resultGroupsDown <- read_tsv("pathways2experiments_significant_Down_0.05.txt") %>% rename(X1 = "name")
resultGroupsUp <- read_tsv("pathways2experiments_significant_Positive_0.05.txt");rename(resultGroupsUp,replace = c("X1" = "name"))
#resultGroupsUp %<>% filter( grepl("GO:", name)) %>% mutate(goID = gsub(".*%","" ,name))
resultGroupsDown <- read_tsv("pathways2experiments_significant_Negative_0.05.txt"); rename(resultGroupsDown,replace = c("X1" = "name"))
#resultGroupsDown %<>% filter( grepl("GO:", name)) %>% mutate(goID = gsub(".*%","" ,name))
(resultGroups <- bind_rows(resultGroupsUp, resultGroupsDown))
#resultGroups %<>% filter( grepl("GO:", name)) %>% mutate(goID = gsub(".*%","" ,name))
(resultGroups %<>% filter( grepl("GO:", X1)) %>% mutate(goID = gsub(".*%","" ,X1)))
#write.csv(resultGroups,file = "AgeCT_SignificantPathways_0.05.csv")

#resultGroups <- read.csv("CT_SignificantPathways_0.05.csv",stringsAsFactors = F)

GOAncestorMap <- as.list(GOBPANCESTOR)
#GOAncestorMap <- as.list(GOMFANCESTOR)
#GOAncestorMap <- as.list(GOCCANCESTOR)

# Remove GO IDs that do not have any ancestor (there are none!)
# Match GOID from the results and filter it
# Sort in decreasing order
# Make a dataframe of ancestor count for all the GOIDs in our result
# Rename the column as "occuranceInTable"
# Add one more column of GOID (why?)
# Filter out "all"
GOAncestorMap <- GOAncestorMap[!is.na(GOAncestorMap)] 
goAncestorCounts <- table(unlist(GOAncestorMap[resultGroups$goID])) 
goAncestorCounts <- sort(goAncestorCounts, decreasing = T ) 
df <- data.frame(matrix(unlist(goAncestorCounts), nrow=length(goAncestorCounts), byrow=T), row.names = names(goAncestorCounts))
head(df)
colnames(df)[1] <- "occuranceInTable"
df$GOID <- names(goAncestorCounts)
head(df)
df <- tbl_df(df) %>% filter(GOID != "all") 
head(df)
df <- mutate(rowwise(df), ancestorCount = length(unlist(GOAncestorMap[GOID])))
head(df)
df <- mutate(rowwise(df), name = Term(GOID)) %>% arrange(ancestorCount, desc(occuranceInTable))
head(df)
df$percentOfGOgroups <- df$occuranceInTable/nrow(resultGroups) * 100
head(df)
t(GOAncestorMap[resultGroups$goID])

list.files(pattern = "*.csv")
selectedGroups <- read_tsv("SelectedGroupsBP.tsv")
#selectedGroups <- read_tsv("SelectedGroupsCC.tsv")
#selectedGroups <- read_tsv("SelectedGroupsMF.tsv")

##################################################
resultGroups <- as.data.frame(resultGroups)
head(rownames(resultGroups))
# For each selectedGOID that is in selectedGroups$GOID print the GOID and the associated GO terms
for(selectedGOID in selectedGroups$GOID) { #for each targetted GO group
  print(selectedGOID)
  print(Term(selectedGOID))
  }

# For each rowNumber that is in resultGroups (and given by commands rownames(resultGroups)) print the GOID and the associated GO terms

for(rowNumber in rownames(resultGroups)) { #for each row in the result table
    print(goID <- resultGroups[rowNumber, "goID"])
    print(rowName <- Term(selectedGOID))
}
    if(is.na(rowName)) { rowName <- selectedGOID }
    resultGroups[rowNumber, rowName]= selectedGOID %in% unlist(GOAncestorMap[goID]) #test if it is one of the ancestors/parents
#######################  




resultGroups <- as.data.frame(resultGroups)
head(rownames(resultGroups))
for(selectedGOID in selectedGroups$GOID) { #for each targetted GO group
  print(selectedGOID)
  print(Term(selectedGOID))
  for(rowNumber in rownames(resultGroups)) { #for each row in the result table
    goID <- resultGroups[rowNumber, "goID"]
    rowName <- Term(selectedGOID)
    if(is.na(rowName)) { rowName <- selectedGOID }
    resultGroups[rowNumber, rowName]= selectedGOID %in% unlist(GOAncestorMap[goID]) #test if it is one of the ancestors/parents
  }
}

write.csv(resultGroups,"Sandbox_TruthTable_BP.csv")

expVars <- c("PYC","PV","SST","VIP")
allResults <- NULL
for(expVar in expVars) {
  resultGroupsSubSet <- resultGroups
  resultGroupsSubSet <- resultGroupsSubSet[,c(setdiff(colnames(resultGroupsSubSet), expVars))]
  meltResult <- tbl_df(melt(id.vars = c("X1", "goID"), resultGroupsSubSet))
  meltResultNumberTrue <- meltResult %>% group_by(variable) %>% dplyr::summarise(n=sum(value))
  for (direction in c(1,-1)) {
    print(direction)
    #filter so that the desiered sign is the same as the effect
    resultGroupsSubSet <- resultGroups[sign(direction) ==  sign(resultGroups[,expVar]),]
    resultGroupsSubSet <- resultGroupsSubSet[,c(setdiff(colnames(resultGroupsSubSet), expVars))]
    meltResult <- tbl_df(melt(id.vars = c("X1", "goID"), resultGroupsSubSet))
    meltResult <- meltResult %>% group_by(variable) %>% dplyr::summarise(directionTrue=sum(value))
    meltResult <- inner_join(meltResultNumberTrue,meltResult) %>% mutate(proportion = directionTrue/n) %>% dplyr::select(variable, proportion)
    colnames(meltResult) <- c("GO Group", paste0(expVar, ifelse(direction==1, ".up", ".down")))
    if(is.null(allResults) || nrow(allResults)==0) {
      allResults <- meltResult
    } else {
      if (nrow(meltResult) != 0) {
        allResults <- inner_join(allResults, meltResult)
      }
    }
  }
}



write.csv(allResults,"Samdbox_SummaryBP.csv")
