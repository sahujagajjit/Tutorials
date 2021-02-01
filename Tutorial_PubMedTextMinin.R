## Date: 1 February 2021 ######################################################
## Author: Jagajjit Sahu ######################################################
## Mining PubMed Texts ########################################################
###############################################################################
######################################START####################################
###############################################################################

# Load library RISmed
library(RISmed)
# Set working directory
setwd("E:/Mine/Presentations/Anthony's_2021/30Jan21")
# Create a folder and set path to it
dir.create("PubMedMinin")
setwd("PubMedMinin")
# A way of doing it is to download all PMIDs and then use package RISmed to download all title and abstract
my_query <- EUtilsSummary("Jagajjit Sahu[AU]",retmax=50)
summary(my_query)
# get actual data from PubMed
records<- EUtilsGet(my_query)
# Extract PMIDs
PMID(records)
# Extract Year of publication
YearPubmed(records)
# Create a data frame for PMID, Title, Abstract and Year
pubmedDat <- data.frame(PMID=PMID(records), Title=ArticleTitle(records), Abstract=AbstractText(records), Year=YearPubmed(records))
# Save the dataframe
write.table(pubmedDat, file = "pmidzTitleAbstractYear.txt", sep = "\t" , row.names = F)
# Plot for Year wise publication
## Find count of publication for each year
yearWiseCount <- table(pubmedDat$Year)
## Check the class of the object yearWiseCount
class(yearWiseCount)
## Convert to dataframe
yearWiseCount <- as.data.frame(yearWiseCount)
## Change the columnnames
### Check the column names
colnames(yearWiseCount)
### Rename the column names
colnames(yearWiseCount) <- c("YearOfPublication", "NoOfPublications")
## Plot
barplot(yearWiseCount$NoOfPublications)
library(ggpubr)
ggbarplot(yearWiseCount, x = "YearOfPublication", y = "NoOfPublications",
          fill = "YearOfPublication", color = "YearOfPublication", palette = "jco")

# Read the saved file to create a word cloud from the frequent words in the Title and Abstract
pmidzTitleAbstractYear <- read.table("pmidzTitleAbstractYear.txt", header = T, sep = "\t", stringsAsFactors = F, encoding = "UTF-8")

library("tm")
library("wordcloud")

# Text processing (Removing stop words, etc.) and creating wordcloud
stopWords <- stopwords("en")
# Title and Abstract merged for wordcloud
titlAbst <- rbind(pmidzTitleAbstractYear$Title, pmidzTitleAbstractYear$Abstract)
# Creating a corpus from the text selected
titlAbstDocs <- Corpus(VectorSource(titlAbst))
# Creating a content transfermer function to preprocess the text
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
# Using the replacement function above with tm_map function to get rid of the noises as symbols
titlAbstDocs <- tm_map(titlAbstDocs, toSpace, "/")
titlAbstDocs <- tm_map(titlAbstDocs, toSpace, "@")
titlAbstDocs <- tm_map(titlAbstDocs, toSpace, "\\|")
# Convert the text to lower case
titlAbstDocs <- tm_map(titlAbstDocs, content_transformer(tolower))
# Remove numbers
titlAbstDocs <- tm_map(titlAbstDocs, removeNumbers)
# Remove english common stopwords
titlAbstDocs <- tm_map(titlAbstDocs, removeWords, stopWords)
# Remove your own stop word
## specify your stopwords as a character vector
titlAbstDocs <- tm_map(titlAbstDocs, removeWords, c("blabla1", "blabla2")) 
# Remove punctuations
titlAbstDocs <- tm_map(titlAbstDocs, removePunctuation)
# Eliminate extra white spaces
titlAbstDocs <- tm_map(titlAbstDocs, stripWhitespace)
# Text stemming
# titlDocs <- tm_map(titlDocs, stemDocument)
titlAbstDtm <- TermDocumentMatrix(titlAbstDocs)
titlAbstDtmMat <- as.matrix(titlAbstDtm)
titlAbstDtmMatSortd <- sort(rowSums(titlAbstDtmMat),decreasing=TRUE)
titlAbstDtmMatSortd.df <- data.frame(word = names(titlAbstDtmMatSortd),freq=titlAbstDtmMatSortd)
write.csv(titlAbstDtmMatSortd.df, file = "titlAbstDtmMatSortd.df.csv", row.names = F)
# Save the wordcloud as a png image
png("titlAbstWordcloud.png", height=1000, width=1000)
wordcloud(words = titlAbstDtmMatSortd.df$word, freq = titlAbstDtmMatSortd.df$freq, min.freq = 1,
          max.words=500, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
dev.off()
# Remove generic words and plot wordcloud again
genericWords <- c("genes", "sequence", "gene")
titlAbst_genericWordRmvd.df <- titlAbstDtmMatSortd.df[(titlAbstDtmMatSortd.df$word != genericWords),]
# Save the wordcloud as a png image
png("titlAbstWordcloudGenrWOrdsRmvd.png", height=1000, width=1000)
wordcloud(words = titlAbst_genericWordRmvd.df$word, freq = titlAbst_genericWordRmvd.df$freq, min.freq = 1,
          max.words=500, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
dev.off()

###############################################################################
######################################END######################################
###############################################################################