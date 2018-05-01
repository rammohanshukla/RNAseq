#pval_thresh <- 0.001
pval_thresh <- 0.05
fdr_thresh <- 0.25
NES_thresh <- 0
min_experiments <- 1

similarity_cutoff <- "0.5"
similarity_metric <- "OVERLAP"

#set working directory
#setwd("./RAM_Etiennelab_celltypes/");

gsea_enrichments <- read.table("Gsea_reports_allcelltypes.txt", header = TRUE, sep = "\t",as.is = TRUE, quote="\"")

colnames(gsea_enrichments) <- c("Experiment","Name","GS","GS DETAILS","SIZE","ES", "NES",
                                "NOM p-val","FDR q-val","FWER p-val","RANK AT MAX","LEADING EDGE");


#filter gsea enrichment by thresholds
#only include NES scores that are > 0 (for the proteomics data the under-enriched is not fitting for this analysis)
gsea_enrichments_filtered <- gsea_enrichments[which(gsea_enrichments[,'NOM p-val']<=pval_thresh & 
                                                    gsea_enrichments[,'FDR q-val']<=fdr_thresh & 
                                                    gsea_enrichments[,'NES']<NES_thresh),]


#Save the two files
#write.csv(gsea_enrichments,file = "gsea_enrichments.csv")
write.csv(gsea_enrichments_filtered,file = "gsea_enrichments_filtered_Down.csv")
#row_names - get the unique set of pathways that is contained in the collated data.  column 2 indicates the pathway name
row_names <- unique(gsea_enrichments_filtered[,2])
length(row_names)
#column_names - get the unique set of experiments contains in the collated data.  column 1 indicates the experiment type
column_names <- unique(gsea_enrichments_filtered[,1])
length(column_names)
#colnames(pathways2experiments_significant) <- c("HypoH","HypoA","HypoB","HypoC","HypoF","HypoG")

#create a matrix which will store the cell type profiles for all genesets in the thresholded set
pathways2experiments_significant <- matrix(nrow=length(row_names), ncol=length(column_names),dimnames=list(row_names, column_names))
pathways2experiments_all <- matrix(nrow=length(row_names), ncol=length(column_names),dimnames=list(row_names, column_names))

for (i in 1:length(row_names)){
  for (j in 1:length(column_names)){
      #only add the NES value if this pathway is significant for this experiment
    if(length(which(gsea_enrichments_filtered[,1] == column_names[j] & 
                    gsea_enrichments_filtered[,2] == row_names[i])) > 0 ){
        pathways2experiments_significant[i,j] = gsea_enrichments_filtered[which(gsea_enrichments_filtered[,1] == column_names[j] & 
                                                               gsea_enrichments_filtered[,2] == row_names[i]), 'NES']   
    }
      #only add the NES value if it exists in the enrichments (irrespective of significance)
      if(length(which(gsea_enrichments[,1] == column_names[j] & 
                    gsea_enrichments[,2] == row_names[i])) > 0 ){
      pathways2experiments_all[i,j] = gsea_enrichments[which(gsea_enrichments[,1] == column_names[j] & 
                                                               gsea_enrichments[,2] == row_names[i]), 'NES']
      }
  }
}

# only include the pathways that are significant in at least X experiments

pathways2experiments_significant <- pathways2experiments_significant[
    which(apply(pathways2experiments_significant,1, function(x){length(which(x!=0))}) >= min_experiments),]
pathways2experiments_all <- pathways2experiments_all[
    which(apply(pathways2experiments_significant,1, function(x){length(which(x!=0))}) >= min_experiments),]

# create a fake gsea summary enrichment file.  Can't use generic without loading in the gmt file and computing genes
# that belong to each set.  In order to use gsea option enrichment file needs to be in the same format as gsea.

#create a collapsed enrichment file
# NES and ES are just the summed values
# p-value and q-value is the average
collapsedenr_column_names = c("Name","GS","GS DETAILS","SIZE","ES", "NES","NOM p-val","FDR q-val","FWER p-val","RANK AT MAX","LEADING EDGE");

#limit to a subset
row_names_subset <- rownames(pathways2experiments_significant)[
    which(apply(pathways2experiments_significant,1, function(x){length(which(x!=0))}) >= min_experiments)]

collapsed_enrichments<- matrix(nrow=length(row_names_subset), ncol=length(collapsedenr_column_names),
                               dimnames=list(row_names_subset, collapsedenr_column_names))
#go through the genesets
for (i in 1:length(row_names_subset)){
  #get all the genesets from the filtered set
  indices <- which(gsea_enrichments_filtered[,2] == row_names_subset[i])
  subset <- gsea_enrichments_filtered[indices,];
  #grab the first name,gs,gs details, size - they are all the same 
  collapsed_enrichments[i,1] <- subset[1,2]
  collapsed_enrichments[i,2] <- subset[1,3]
  collapsed_enrichments[i,3] <- subset[1,4]
  collapsed_enrichments[i,4] <- subset[1,5]
  #get the summed ES score
  collapsed_enrichments[i,5] <- sum(subset[,6])
  #get the summed NES score
  collapsed_enrichments[i,6] <- sum(subset[,7])
  #get the average pvalue
  collapsed_enrichments[i,7] <- 0.05
  #get the average fdr
  collapsed_enrichments[i,8] <- 0.05
  #get the average FWER
  collapsed_enrichments[i,9] <- mean(subset[,10])
  #get the first rank at max
  collapsed_enrichments[i,10] <- max(subset[,11])
  collapsed_enrichments[i,11] <- subset[1,12]

}


collapsed_enrichments[,1] <- trimws(collapsed_enrichments[,1])
collapsed_enrichments[,2] <- trimws(collapsed_enrichments[,2])
rownames(collapsed_enrichments) <- trimws(rownames(collapsed_enrichments))

#output the pathways2experiments_all
#
# The pathways2experiments file is a matrix of pathways to experiments where each value in the matrix is a significant NES value.
# This table can be used to calculate which genesets pass the minimum expereiment threshold but should not be used 
# as an expression file for vista clara plugin as it is missing NES values for pathways and experiments that were not significant

pathways2experiments_all[is.na(pathways2experiments_all)] <- 0
#we are not interested in the sets that are significant in rest of the tissues.
pathways2experiments_all[which(pathways2experiments_all < 0)] <- 0
rownames(pathways2experiments_all) <- trimws(rownames(pathways2experiments_all))
output_pathways2experiments_all <- cbind(rownames(pathways2experiments_all), pathways2experiments_all)

output_pathways2experiments_all <- output_pathways2experiments_all[grep(output_pathways2experiments_all[,1],
                                              pattern = "TCR SIGNALING IN NA",
                                              invert=TRUE),]
output_pathways2experiments_all<- output_pathways2experiments_all[grep(output_pathways2experiments_all[,1],
                                              pattern = "LOSS OF PROTEINS REQUIRED FOR INTERPHASE MICROTUBULE",
                                              invert=TRUE),]
output_pathways2experiments_all<- output_pathways2experiments_all[grep(output_pathways2experiments_all[,1],
                                              pattern = "DOWNSTREAM SIGNALING IN NA",
                                              invert=TRUE),]
write.table(output_pathways2experiments_all, file=paste0("pathways2experiments_all_Down",pval_thresh,".txt",sep=""), sep="\t", 
            row.names=FALSE, col.names=TRUE,quote=FALSE)


#output the pathways2experiments_significant as well
#

pathways2experiments_significant[is.na(pathways2experiments_significant)] <- 0
rownames(pathways2experiments_significant) <- trimws(rownames(pathways2experiments_significant))
output_pathways2experiments_significant <- cbind(rownames(pathways2experiments_significant), pathways2experiments_significant)

output_pathways2experiments_significant <- output_pathways2experiments_significant[grep(output_pathways2experiments_significant[,1],
                                              pattern = "TCR SIGNALING IN NA",
                                              invert=TRUE),]
output_pathways2experiments_significant<- output_pathways2experiments_significant[grep(output_pathways2experiments_significant[,1],
                                              pattern = "LOSS OF PROTEINS REQUIRED FOR INTERPHASE MICROTUBULE",
                                              invert=TRUE),]
output_pathways2experiments_significant <- output_pathways2experiments_significant[grep(output_pathways2experiments_significant[,1],
                                              pattern = "DOWNSTREAM SIGNALING IN NA",
                                              invert=TRUE),]

write.table(output_pathways2experiments_significant, file=paste0("pathways2experiments_significant_",pval_thresh,".txt",sep=""), sep="\t", 
            row.names=FALSE, col.names=TRUE,quote=FALSE)


library(pheatmap)
library("RColorBrewer")
col.pal <- rev(brewer.pal(3,"Reds"))

new_names <- lapply(rownames(pathways2experiments_significant),function(x){unlist(strsplit(x,"%"))[1]} )

dim(pathways2experiments_significant)
rownames(pathways2experiments_significant) <- new_names
tree_results <- pheatmap(mat=pathways2experiments_significant,show_rownames = FALSE,color = col.pal,)

############

colnames(pathways2experiments_significant)
rownames(pathways2experiments_significant)[which(pathways2experiments_significant[, "OY_SST"]<0 &
      pathways2experiments_significant[, "OY_PYC"] == 0 &
      pathways2experiments_significant[, "OY_VIP"] == 0 &
      pathways2experiments_significant[, "OY_PV"] == 0)]

cutree(tree_results$tree_row,k=4)

#to save on computation you can reload previous computed collapsed enrichments
#enrichment_results_file_name <- paste("mastermap_enrichments_fdr", fdr_thresh, "_minexp", min_experiments,".txt", sep="")
#collapsed_enrichments <- read.table(enrichment_results_file_name, , header = TRUE, sep = "\t",as.is = TRUE, quote="\"")

Sys.getlocale()
Sys.setlocale(locale="C")

collapsed_enrichments <- collapsed_enrichments[grep(collapsed_enrichments[,1],
                                              pattern = "TCR SIGNALING IN NA",
                                              invert=TRUE),]
collapsed_enrichments <- collapsed_enrichments[grep(collapsed_enrichments[,1],
                                              pattern = "LOSS OF PROTEINS REQUIRED FOR INTERPHASE MICROTUBULE",
                                              invert=TRUE),]
collapsed_enrichments <- collapsed_enrichments[grep(collapsed_enrichments[,1],
                                              pattern = "DOWNSTREAM SIGNALING IN NA",
                                              invert=TRUE),]

enrichment_results_file_name <- paste("mastermap_enrichments_pval", pval_thresh, "_minexp", min_experiments,".txt", sep="")
write.table( collapsed_enrichments, file=enrichment_results_file_name, sep="\t", row.names=FALSE, col.names=TRUE,quote=FALSE)

library(RJSONIO)
library(httr)

port.number=1234
base.url = paste("http://localhost:",toString(port.number),"/v1", sep="")

print(base.url)

version.url = paste(base.url,"version", sep="/")
cytoscape.version = GET(version.url)
cy.version = fromJSON(rawToChar(cytoscape.version$content))
print(cy.version)

enrichmentmap.url <- paste(base.url, "commands", "enrichmentmap", "build", sep="/")
#mac file paths
gmt_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/GSEA_results/Human_AllPathwaysGOBP_noiea_May_14_2013_symbol.gmt"
path_to_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/GSEA_results/"
exp_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/Human_Proteome_Map_spectral_count_gene_tissue.txt"

gmt_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/Mastermap_example/notebooks/RAM_Etiennelab_celltypes/Mouse_GOBP_AllPathways_no_GO_iea_September_01_2016_symbol.gmt"
path_to_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/Mastermap_example/notebooks/RAM_Etiennelab_celltypes/"
exp_file="/Users/risserlin/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/Mastermap_example/notebooks/RAM_Etiennelab_celltypes/oldyoung_expression.txt"
#wi

#windows file paths
#gmt_file="C:/Users/zaphod/Ruth_dropbox/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/GSEA_results/Human_AllPathwaysGOBP_noiea_May_14_2013_symbol.gmt"
#exp_file="C:/Users/zaphod/Ruth_dropbox/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/Human_Proteome_Map_spectral_count_gene_tissue.gct"
#path_to_file="C:/Users/zaphod/Ruth_dropbox/Dropbox (Bader Lab)/Ruth Isserlin's files/Enrichment_Analyses/Mastermap/notebooks/One_expression_file_example/GSEA_results/"

enr_file = paste(path_to_file,enrichment_results_file_name,sep="")

em_params <- list(analysisType = "gsea",enrichmentsDataset1 = enr_file,pvalue="1.0",qvalue="0.05",
                  expressionDataset1 = exp_file, gmtFile = gmt_file,
                  similaritycutoff=similarity_cutoff,coeffecients=similarity_metric)

response <- GET(url=enrichmentmap.url, query=em_params)

content(response, "text",encoding="ISO-8859-1")
