# Loading the original count data
seqData <- read.table('GSE60450_Lactation-GenewiseCounts.txt', header = T, sep = '\t')


# Loading the experimental design file. Each row of expDesign describes the column meaning in seqData. The expDesign is useful for DEseq analysis.
expDesign <- read.table('SampleInfo_Corrected.txt', header = T, sep = '\t')
expDesign

## loading DEseq2 
library(DESeq2)
library(dbplyr)
## This function is designed for differential expression calling in the current study
## The function (deSeqFun) input is: 
## 1. reads count file (row name should be the geneID)
## 2. vector corresponding the design information (e.g. treatment vs. control)


deSeqFun <- function(dfInputCounts, 
                     con = c(rep('c',4), rep('t',2)))
{
  #Calculate the Counts Per Million measure
  myCPM     <- apply(dfInputCounts, 2, function(x){x/sum(x) * 10^6})
  #Identify genes with at least 0.5 cpm in at least 2 samples
  thresh    <- myCPM > 0.5
  keep      <- rowSums(thresh) >= 2
  dfInput   <- dfInputCounts[keep, ]
  
  table2    <- data.frame(name = colnames(dfInput),
                          condition = as.factor(con))
  dds       <- DESeqDataSetFromMatrix(dfInput, 
                                      colData=table2, design= ~ condition)
  dds       <- DESeq(dds)
  norCounts <- counts(dds, normalized=TRUE)
  res       <- results(dds)
  idNew     <- match(rownames(dfInputCounts), rownames(norCounts))
  
  norCountsReturn <- norCounts[idNew, ]
  resReturn       <- res[idNew, ]
  rownames(norCountsReturn) <- rownames(dfInputCounts)
  rownames(resReturn)       <- rownames(dfInputCounts)
  return(list(nC = norCountsReturn, res=resReturn))
}
rownames(seqData) <- seqData$EntrezGeneID
## Reads count extracting for virgin-like state
df_virgin  <- seqData[,grep('virgin', expDesign$Status)  +2 ]
## Reads count extracting for pregnancy state
df_pregnant<- seqData[,grep('pregnant', expDesign$Status)+2 ]
## Reads count extracting for lactation state
df_lactate <- seqData[,grep('lactate', expDesign$Status) +2 ]
## All the matrix here have the same dimension, and the beginning two columns are 'basal,' while the following two columns are "luminal"

## Differential expression calling for virgin-like state
deseq_virgin  <- deSeqFun(df_virgin  ,c('b','b', 'v', 'v'))
## Differential expression calling for pregnancy state
deseq_pregn   <- deSeqFun(df_pregnant,c('b','b', 'v', 'v'))
##Differential expression calling for lactation state
deseq_lactate  <- deSeqFun(df_lactate ,c('b','b', 'v', 'v'))
## Now, the differential expressed genes between basal and luminal for three stages are returned
## The return elements inside each variable include
## 1. [nC], normalized DEseq counts
## 2. [res],statistics value from DEseq
## e.g. Using deseq_lactate$res or deseq_lactate$nC to check 

## Extracting differential expressed gene with padj < 0.01
id_sig_virgin <- which(deseq_virgin$res$padj < 0.01) 
id_sig_pregn  <- which(deseq_pregn$res$padj  < 0.01)
id_sig_lactate <- which(deseq_lactate$res$padj < 0.01)
##Now, only the rowID of differential expressed genes are extracted
# save deseq outputs as csvs
write.csv(deseq_lactate,"deseq_lactate_file.csv")
write.csv(deseq_virgin,"deseq_virgin_file.csv")
write.csv(deseq_pregn,"deseq_preg_file.csv")

## 

virgin2 <- read.csv("deseq_virgin_file.csv")
preg2 <- read.csv("deseq_preg_file.csv")
lactate2 <- read.csv("deseq_lactate_file.csv")


# The first column is emtpy and R assings an X to it we want to rename that 1st  column to "Gene_ID" to be able to merge it later
colnames(virgin2)[1] <- "Gene_ID" 
colnames(lactate2)[1] <- "Gene_ID"
colnames(preg2)[1] <- "Gene_ID"


#filter out values with p adjusted under 0.01
virgin3 <-subset(virgin2, res.padj < 0.01)
#filter out upregulated genes
upregulated_virgin <-subset(virgin3, res.log2FoldChange > 1.0)
#filter out downregulated genes
downregulated_virgin <-subset(virgin3, res.log2FoldChange < -1.0)


#filter out values with p adjusted under 0.01
preg3 <-subset(preg2, res.padj < 0.01)
#filter out upregulated genes
upregulated_preg <-subset(preg3, res.log2FoldChange > 1.0)
#filter out downregulated genes
downregulated_preg <-subset(preg3, res.log2FoldChange < -1.0)

#filter out values with p adjusted under 0.01
lactate3 <-subset(lactate2, res.padj < 0.01)
#filter out upregulated genes
upregulated_lactate <-subset(lactate3, res.log2FoldChange > 1.0)
#filter out downregulated genes
downregulated_lactate <-subset(lactate3, res.log2FoldChange < -1.0)

## differential expressed genes at all stages 
library(gplots)
library(VennDiagram)
library(ggVennDiagram)
# Venn diagram for upregulated genes

list_virgin=as.list(upregulated_virgin$Gene_ID)
list_preg=as.list(upregulated_preg$Gene_ID)
list_lactate=as.list(upregulated_lactate$Gene_ID)

x <- list( 
  A = list_virgin,
  B = list_preg,
  C = list_lactate)

venn(x)


#Venn diagram for downregulated genes

list_virgin_down=as.list(downregulated_virgin$Gene_ID)
list_preg_down=as.list(downregulated_preg$Gene_ID)
list_lactate_down=as.list(downregulated_lactate$Gene_ID)

y <- list( 
  A = list_virgin_down,
  B = list_preg_down,
  C = list_lactate_down)

venn(y)

list 
## Bar plot drawing, showing the number of differential expressed genes
## May considering give a Venn diagram for overlapped differential expressed genes at different stage?
par(mar = c(3,8,1,1))
barplot(c(length(id_sig_virgin),
          length(id_sig_pregn),
          length(id_sig_lactate)), col = 'cornflowerblue',
        names = c('virgin','pregn', 'lactate'), las = 1,
        ylab = '#differential expressed gene', ylim = c(0,10000))

## Normalizing the data by TPM
tpm <- apply(seqData[,-c(1,2)], 2, function(x){
  x/seqData$Length*10^9/sum(x)})
rownames(tpm) <- seqData$EntrezGeneID

## Extracting the rowID of differential expressed gene 
sigRowID  <- unique(c(id_sig_virgin, id_sig_pregn, id_sig_lactate))
tpm_sigDif<- tpm[sigRowID, ]

## Performing PCA
tpm_sigDifPCA   <- prcomp(tpm_sigDif,scale=F )
## Assigning different color to the samples
expDesign$color <- c('red','red',#'basal & virgin'
                     'forestgreen','forestgreen',#'basal & pregnant'
                     'gold','gold', #'basal & lactate'
                     'cornflowerblue','cornflowerblue', #'luminal & virgin'
                     'purple','purple', #'luminal & pregnant'
                     'gray','gray') #'luminal & lactate'
## PCA plot
plot(tpm_sigDifPCA$rotation[,1],
     tpm_sigDifPCA$rotation[,2], 
     bg = expDesign$color,pch = 21,cex = 2,
     las = 1,xlab = 'PCA_direction1',
     ylab = 'PCA_direction2')
legend('bottomright', 
       legend = unique(paste(expDesign$CellType,
                             expDesign$Status, sep ='_')),
       col= c('red','forestgreen','gold',
              'cornflowerblue', 'purple', 'gray'),
       pch = rep(19,6),cex=rep(1,6),
       bty = 'n')

## Generating a data frame,using color indicating up or down regulated genes
res_plot      <- data.frame( deseq_virgin$res)
res_plot$col  <- 'gray40'
res_plot$col[which(res_plot$log2FoldChange > 1 & res_plot$padj < 0.01)] <- 'red'
res_plot$col[which(res_plot$log2FoldChange < -1 & res_plot$padj < 0.01)] <-'cornflowerblue'

plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)', 
     las = 1 
)

res_plot      <- data.frame( deseq_pregn$res)
res_plot$col  <- 'gray40'
res_plot$col[which(res_plot$log2FoldChange > 1 & res_plot$padj < 0.01)] <- 'red'
res_plot$col[which(res_plot$log2FoldChange < -1 & res_plot$padj < 0.01)] <-'cornflowerblue'

plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)', 
     las = 1 
)

res_plot      <- data.frame( deseq_lactate$res)
res_plot$col  <- 'gray40'
res_plot$col[which(res_plot$log2FoldChange > 1 & res_plot$padj < 0.01)] <- 'red'
res_plot$col[which(res_plot$log2FoldChange < -1 & res_plot$padj < 0.01)] <-'cornflowerblue'

plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)', 
     las = 1 
)


## Loading the KEGG annotation file
geneID_KeggPathway <- read.table('kegg.pathway.gene.txt', header = T, sep = '\t')

## Building the funciton for KEGG enrichment analysis
## The function (deSeqFun) input is: 
## 1. Gene ID for test
## 2. All gene ID in background
## 3. Matrix of KEGG pathway annotation. The first column is geneID, while the second column is corresponding KEGG annotation
enrichment_analysis <- function(geneID_query,
                                geneID_background,
                                geneAnnotation) 
{
  geneQueryInf      <- geneAnnotation[geneAnnotation[,1] %in% geneID_query,]
  geneBackgroundInf <- geneAnnotation[geneAnnotation[,1] %in% geneID_background,]
  queryGOterm <- unique(geneQueryInf[,2])
  
  goResult <- c()
  aa <- sapply(queryGOterm, function(id) 
  {
    numQuery        <- length( which(geneQueryInf[,2] == id) )
    numBackground   <- length( which(geneBackgroundInf[,2] == id))
    
    numQuery_no     <- length(geneID_query) - numQuery
    numBackground_no<- length(geneID_background) - numBackground
    
    #print(c(numBackground, numBackground_no))
    fishTest <- fisher.test(rbind( c(numQuery, numQuery_no),
                                   c(numBackground, numBackground_no) ),
                            alternative = 'greater')
    infReturn <- c(numQuery,
                   numQuery_no,
                   numBackground,
                   numBackground_no, fishTest$p.value)
    goResult <<- rbind(goResult, infReturn)
  })
  rownames(goResult) <- queryGOterm
  colnames(goResult) <- c('#QueryWithKEGGterm',
                          '#QueryWithoutKEGGterm',
                          '#BackgroundWithKEGGterm',
                          '#BackgroundWithoutKEGGterm', 'pvalue')
  goResult <- data.frame(goResult)
  goResult$padj <- p.adjust(goResult$pvalue, method = 'fdr')
  return(goResult)
}
## Extracting differential expressed genes during lactation
dif_lactateGene <- rownames(deseq_lactate$res)[which(deseq_lactate$res$padj < 0.01)]
## Performing enrichment analysis using function enrichment_analysis
keggEnrichment_lactate <- enrichment_analysis(
  dif_lactateGene[dif_lactateGene %in% geneID_KeggPathway$geneID],   ## query geneID in whole KEGG list
  unique(geneID_KeggPathway$geneID),     ## geneID in whole KEGG list
  geneID_KeggPathway)                    ## KEGG list, two columns 


## Picking up the KEGG term with padj < 0.1 for state lactation
keggSig_lactate   <- keggEnrichment_lactate[keggEnrichment_lactate$padj < 0.1,]
## Calculating the expection of gene counts. Based on the fraction of the genes with specific KEGG term, multiple the total number of the genes in query data sets (e.g. lactation state)
keggSig_lactate$expection <- keggSig_lactate$X.BackgroundWithKEGGterm/(keggSig_lactate$X.BackgroundWithKEGGterm +
                                                                         keggSig_lactate$X.BackgroundWithoutKEGGterm) * (keggSig_lactate$X.QueryWithKEGGterm + 
                                                                                                                           keggSig_lactate$X.QueryWithoutKEGGterm)

keggSig_lactate <- keggSig_lactate[order(keggSig_lactate$X.QueryWithKEGGterm), ]
keggSigDraw_lactate <- t( keggSig_lactate[,c(1,7)] )

grep('mTOR', rownames(keggEnrichment_lactate), ignore.case = T)

## Bar plot drawing
par(mar = c(4,20,1,1))
barplot(keggSigDraw_lactate, horiz = T, las = 1)


## Loading the KEGG annotation file
geneID_KeggPathway <- read.table('kegg.pathway.gene.txt', header = T, sep = '\t')

## Building the funciton for KEGG enrichment analysis
## The function (deSeqFun) input is: 
## 1. Gene ID for test
## 2. All gene ID in background
## 3. Matrix of KEGG pathway annotation. The first column is geneID, while the second column is corresponding KEGG annotation
enrichment_analysis <- function(geneID_query,
                                geneID_background,
                                geneAnnotation) 
{
  geneQueryInf      <- geneAnnotation[geneAnnotation[,1] %in% geneID_query,]
  geneBackgroundInf <- geneAnnotation[geneAnnotation[,1] %in% geneID_background,]
  queryGOterm <- unique(geneQueryInf[,2])
  
  goResult <- c()
  aa <- sapply(queryGOterm, function(id) 
  {
    numQuery        <- length( which(geneQueryInf[,2] == id) )
    numBackground   <- length( which(geneBackgroundInf[,2] == id))
    
    numQuery_no     <- length(geneID_query) - numQuery
    numBackground_no<- length(geneID_background) - numBackground
    
    #print(c(numBackground, numBackground_no))
    fishTest <- fisher.test(rbind( c(numQuery, numQuery_no),
                                   c(numBackground, numBackground_no) ),
                            alternative = 'greater')
    infReturn <- c(numQuery,
                   numQuery_no,
                   numBackground,
                   numBackground_no, fishTest$p.value)
    goResult <<- rbind(goResult, infReturn)
  })
  rownames(goResult) <- queryGOterm
  colnames(goResult) <- c('#QueryWithKEGGterm',
                          '#QueryWithoutKEGGterm',
                          '#BackgroundWithKEGGterm',
                          '#BackgroundWithoutKEGGterm', 'pvalue')
  goResult <- data.frame(goResult)
  goResult$padj <- p.adjust(goResult$pvalue, method = 'fdr')
  return(goResult)
}
## Extracting differential expressed genes during virgin
dif_virginGene <- rownames(deseq_virgin$res)[which(deseq_virgin$res$padj < 0.01)]
## Performing enrichment analysis using function enrichment_analysis
keggEnrichment_virgin <- enrichment_analysis(
  dif_virginGene[dif_virginGene %in% geneID_KeggPathway$geneID],   ## query geneID in whole KEGG list
  unique(geneID_KeggPathway$geneID),     ## geneID in whole KEGG list
  geneID_KeggPathway)                    ## KEGG list, two columns 


## Picking up the KEGG term with padj < 0.1 for state virgin
keggSig_virgin   <- keggEnrichment_virgin[keggEnrichment_virgin$padj < 0.1,]
## Calculating the expection of gene counts. Based on the fraction of the genes with specific KEGG term, multiple the total number of the genes in query data sets (e.g. virgin state)
keggSig_virgin$expection <- keggSig_virgin$X.BackgroundWithKEGGterm/(keggSig_virgin$X.BackgroundWithKEGGterm +
                                                                         keggSig_virgin$X.BackgroundWithoutKEGGterm) * (keggSig_virgin$X.QueryWithKEGGterm + 
                                                                                                                           keggSig_virgin$X.QueryWithoutKEGGterm)

keggSig_virgin <- keggSig_virgin[order(keggSig_virgin$X.QueryWithKEGGterm), ]
keggSigDraw_virgin <- t( keggSig_virgin[,c(1,7)] )

grep('mTOR', rownames(keggEnrichment_pregn), ignore.case = T)

## Bar plot drawing
par(mar = c(4,20,1,1))
barplot(keggSigDraw_virgin, horiz = T, las = 1)


## Loading the KEGG annotation file
geneID_KeggPathway <- read.table('kegg.pathway.gene.txt', header = T, sep = '\t')

## Building the funciton for KEGG enrichment analysis
## The function (deSeqFun) input is: 
## 1. Gene ID for test
## 2. All gene ID in background
## 3. Matrix of KEGG pathway annotation. The first column is geneID, while the second column is corresponding KEGG annotation
enrichment_analysis <- function(geneID_query,
                                geneID_background,
                                geneAnnotation) 
{
  geneQueryInf      <- geneAnnotation[geneAnnotation[,1] %in% geneID_query,]
  geneBackgroundInf <- geneAnnotation[geneAnnotation[,1] %in% geneID_background,]
  queryGOterm <- unique(geneQueryInf[,2])
  
  goResult <- c()
  aa <- sapply(queryGOterm, function(id) 
  {
    numQuery        <- length( which(geneQueryInf[,2] == id) )
    numBackground   <- length( which(geneBackgroundInf[,2] == id))
    
    numQuery_no     <- length(geneID_query) - numQuery
    numBackground_no<- length(geneID_background) - numBackground
    
    #print(c(numBackground, numBackground_no))
    fishTest <- fisher.test(rbind( c(numQuery, numQuery_no),
                                   c(numBackground, numBackground_no) ),
                            alternative = 'greater')
    infReturn <- c(numQuery,
                   numQuery_no,
                   numBackground,
                   numBackground_no, fishTest$p.value)
    goResult <<- rbind(goResult, infReturn)
  })
  rownames(goResult) <- queryGOterm
  colnames(goResult) <- c('#QueryWithKEGGterm',
                          '#QueryWithoutKEGGterm',
                          '#BackgroundWithKEGGterm',
                          '#BackgroundWithoutKEGGterm', 'pvalue')
  goResult <- data.frame(goResult)
  goResult$padj <- p.adjust(goResult$pvalue, method = 'fdr')
  return(goResult)
}
## Extracting differential expressed genes during pregn
dif_pregnGene <- rownames(deseq_pregn$res)[which(deseq_pregn$res$padj < 0.01)]
## Performing enrichment analysis using function enrichment_analysis
keggEnrichment_pregn <- enrichment_analysis(
  dif_pregnGene[dif_pregnGene %in% geneID_KeggPathway$geneID],   ## query geneID in whole KEGG list
  unique(geneID_KeggPathway$geneID),     ## geneID in whole KEGG list
  geneID_KeggPathway)                    ## KEGG list, two columns 


## Picking up the KEGG term with padj < 0.1 for state pregn
keggSig_pregn   <- keggEnrichment_pregn[keggEnrichment_pregn$padj < 0.1,]
## Calculating the expection of gene counts. Based on the fraction of the genes with specific KEGG term, multiple the total number of the genes in query data sets (e.g. pregn state)
keggSig_pregn$expection <- keggSig_pregn$X.BackgroundWithKEGGterm/(keggSig_pregn$X.BackgroundWithKEGGterm +
                                                                     keggSig_pregn$X.BackgroundWithoutKEGGterm) * (keggSig_pregn$X.QueryWithKEGGterm + 
                                                                                                                     keggSig_pregn$X.QueryWithoutKEGGterm)

keggSig_pregn <- keggSig_pregn[order(keggSig_pregn$X.QueryWithKEGGterm), ]
keggSigDraw_pregn <- t( keggSig_pregn[,c(1,7)] )

grep('mTOR', rownames(keggEnrichment_virgin), ignore.case = T)
  
## Bar plot drawing
par(mar = c(4,20,1,1))
barplot(keggSigDraw_pregn, horiz = T, las = 1)

print(keggEnrichment_lactate)

