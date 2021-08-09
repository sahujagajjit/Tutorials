# Script for finding miRNA targets, Adding Gene symbols as a column and merging results for number of miRNAs
# Author: JeeT
# Date Created: 03-08-2021
# Last Modified: 09-08-2021

# Define path 
#setwd("E:/Others/")
# Installing packages
#install.packages("BiocManager")
#BiocManager::install("miRNAtap")
#BiocManager::install("miRNAtap.db")
#BiocManager::install("topGO")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")

library(dplyr)
library(miRNAtap)
library(miRNAtap.db)
library(topGO)

# Create a function for more than one miRNA at a time
# Input miRNA IDs; Modify based on the input source such as from file
mirList <- c("miR-10b", "miR-103", "miR-130")
# Define variable for saving all targets
allTargts <- rbind()
for (i in 1:length(mirList)) {
  mir <- mirList[i]
  targts <- getPredictedTargets(mir, species = 'hsa',method = 'geom', min_src = 2)
  if(!is.null(targts)){
    targts <- as.data.frame(targts)
    # Convert the rownames into a column
    targts$EntrezIDs <- rownames(targts)
    # Remove rownames
    rownames(targts) <- NULL
    # Add another column for miRNA
    targts$miRNA <- mir
    allTargts <- bind_rows(allTargts, targts)
  }else{
    print(paste0("No targets were found for: ", mir))
  }
  
}
# Arrange columns
allTargts <- allTargts[, c(ncol(allTargts), (ncol(allTargts)-1), 1:(ncol(allTargts)-2))]
# Save to file
#write.csv(allTargts, file = "allTargts.csv", row.names = F)

# Adding another column as GeneSymbol which should be obtained from the Entrez IDs
library(org.Hs.eg.db)
library(annotate)

# Define the column and get symbols
allTargts$GeneSymbols <- getSYMBOL(allTargts$EntrezIDs, data='org.Hs.eg')
allTargts <- allTargts[, c(1:2, ncol(allTargts), 3:(ncol(allTargts)-1))]
# Save to file
write.csv(allTargts, file = "allTargts.csv", row.names = F)
