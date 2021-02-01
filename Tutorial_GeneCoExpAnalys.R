## Date: 1 February 2021 ######################################################
## Author: Jagajjit Sahu ######################################################
## Gene Expression data Analysis ##############################################
###############################################################################
######################################START####################################
###############################################################################

# set working directory
setwd("E:/Mine/Presentations/Anthony's_2021/30Jan21/GeneExpAnalys")
# Read series matrix data
GSE68646_SeriesMatrix <- read.csv("GSE68646_series_matrixXtrctd.csv", header = T, stringsAsFactors = F)
# Read top table data
GSE68646.top.table <- read.csv("GSE68646.top.table.csv", header = T, stringsAsFactors = F)
# Subset data based on the cutt-off for Fold Change and P-value
## log FC cutt-off 1
GSE68646.top.table_fitrd_FC <- GSE68646.top.table[(GSE68646.top.table$logFC < -1|GSE68646.top.table$logFC > 1),]
## p value cutt-off 0.5
GSE68646.top.table_fitrd_FC_Pval <- GSE68646.top.table_fitrd_FC[(GSE68646.top.table_fitrd_FC$P.Value < 0.5),]
## p value cutt-off 0.1
GSE68646.top.table_fitrd_FC_Pval <- GSE68646.top.table_fitrd_FC[(GSE68646.top.table_fitrd_FC$P.Value < 0.1),]
# Merge gene names and prepare the new series matrix file for the selected ones
## create a subset containing only columns "ID" and "Gene.symbol"
probeGenesymb.df <- GSE68646.top.table_fitrd_FC_Pval[which(colnames(GSE68646.top.table_fitrd_FC_Pval) %in% c("ID", "Gene.symbol"))]
## change the column name ID to ID_REF to merge easily
colnames(probeGenesymb.df)[1] <- "ID_REF"
## merge series matrix and gene symbols dataframe based on the probe ID column
GSE68646_SeriesMatrix_Genesymb <- merge(probeGenesymb.df, GSE68646_SeriesMatrix, by="ID_REF")
# Remove the rows with no Gene symbol
GSE68646_SeriesMatrix_Genesymb_RmvNull <- GSE68646_SeriesMatrix_Genesymb[GSE68646_SeriesMatrix_Genesymb$Gene.symbol != "",]
# Merge rows with same Gene symbol
GSE68646_SeriesMatrix_Final <- aggregate(GSE68646_SeriesMatrix_Genesymb_RmvNull[,-(1:2)], list("Gene.Symbol" = GSE68646_SeriesMatrix_Genesymb_RmvNull[,2]), FUN = "mean")
# Split file into 2 i.e. control and CR
datExprCont <- GSE68646_SeriesMatrix_Final[which(colnames(GSE68646_SeriesMatrix_Final) %in% c("Gene.Symbol", "GSM1678014", "GSM1678015", "GSM1678016", "GSM1678017", "GSM1678018"))]
datExprCR <- GSE68646_SeriesMatrix_Final[which(colnames(GSE68646_SeriesMatrix_Final) %in% c("Gene.Symbol", "GSM1678019", "GSM1678020", "GSM1678021", "GSM1678022", "GSM1678023"))]
# Transpose the dataframe to have genes as column names and the sample names as rownames
geneExprCont <- setNames(data.frame(t(datExprCont[,-1])), datExprCont[,1])
geneExprCR <- setNames(data.frame(t(datExprCR[,-1])), datExprCR[,1])
# Construct co-expression networks for both control and CR
corCont <- cor(geneExprCont)
corCR <- cor(geneExprCR)
## heatmap
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = corCont, col = col, symm = TRUE)
heatmap(x = corCR, col = col, symm = TRUE)
## other analytic plots
### choosing a smaller size data
corCont2 <- corCont[1:30,1:30]
barplot(corCont2)
library(corrplot)
corrplot(corCont2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
library("PerformanceAnalytics")
chart.Correlation(corCont2, histogram=TRUE, pch=19)
# Adjacency matrix
library(igraph)
grphCont <- graph_from_adjacency_matrix(corCont2, mode = "undirected", weighted = T)
grphCR <- graph_from_adjacency_matrix(corCR2, mode = "undirected", weighted = T)
# Cutt off for network
grphCont_cutWgt <- delete.edges(grphCont, which(E(grphCont)$weight < 0.8))
grphCR_cutWgt <- delete.edges(grphCR, which(E(grphCR)$weight < 0.8))
# Intersection of Control and CR graph
grphContCR_cutWgt <- graph.intersection(grphCont_cutWgt, grphCR_cutWgt, keep.all.vertices = FALSE)
# Writing graph to file
write_graph(grphContCR_cutWgt, file = "grphContCR_cutWgt.graphml", format = "graphml")
# Visualization
plot(grphContCR_cutWgt)
library(visNetwork)
visIgraph(grphContCR_cutWgt)

###############################################################################
######################################END######################################
###############################################################################