---
title: "Xenopus figure"
author: "Nick Leigh"
date: "8/3/2022"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
## Load libraries
```{r}
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)


```
## Import data
```{r}
setwd('/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints')
data_dir <- "/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints"

#read in xl data
raw_xl <- Read10X(data.dir = data_dir)

souporcell_calls_s1 <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)


#load metadata from souporcell
souporcell_calls_s1 <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)
```


```{r}
#will simplify these a bit more
souporcell_calls_s1
souporcell_calls_s1[,c(1:2)] -> soup_res_s1

#so now need to see if soup_res when it calls #1 matching up most of the time with one of the CMO calls
dim(soup_res_s1)

#changing column names
names(soup_res_s1) <- c("soup_status", "soup_assign")

```

## Demux with fluors
```{r}
#get flours
fluors <- c('Venus', "mTFP1", "mCherryB51")

#make DF of Fluor counts similar to Cellplex strategy
raw_xl[rownames(raw_xl) %in% fluors,] -> FluorsGenes
as.data.frame(FluorsGenes) -> FluorsDFa
t(FluorsDFa) -> FluorsDF_s1

#make seurat object with the gene expression counts,
#and with souporcell results as the metadata
#so.demux_s1 <- CreateSeuratObject(counts = rna_s1, project = "Round1_spl", assay = 'RNA')
so.demux_s1 <- CreateSeuratObject(counts = raw_xl, meta=souporcell_calls_s1, project = "Xenopus", assay = 'RNA')

#create a new assay for fluorescent counts
fluorlabel_assay <- CreateAssayObject(counts = FluorsGenes)

# add this assay to the previously created Seurat object
so.demux_s1[["fluors"]] <- fluorlabel_assay

#this should output the read depth for each cell. but I need to confirm that this automated
#function is working correctly
so.demux_s1$nCount_RNA -> RNAcounts_s1
as.data.frame(RNAcounts_s1) -> RNAcountsDF_s1
#changing column names in the DF
names(RNAcountsDF_s1) <- c("RNACounts")

#make new column that's the sum of counts for the fluorescent genes for each cell
rowSums(FluorsDF_s1) -> FluorRNAcounts_s1
as.data.frame(FluorRNAcounts_s1) -> FluorRNAcountsDF_s1
#changing column names in the DF
names(FluorRNAcountsDF_s1) <- c("SumFluorCounts")

#run demuxing
#we've found multi works better for us than HTO
so.demux_s1 <- NormalizeData(so.demux_s1, assay = "fluors", normalization.method = "CLR")
so.demux_s1 <- MULTIseqDemux(so.demux_s1, assay = 'fluors', autoThresh = T)
#so.demux$MULTI_ID
length(so.demux_s1$MULTI_ID)
table(so.demux_s1$MULTI_ID)
as.data.frame(so.demux_s1$MULTI_ID) -> multiDF_s1
#changing column names in the DF
names(multiDF_s1) <- c("multi_assign")
so.demux_s1$MULTI_ID-> multi_s1
length(multi_s1)

```

## Compare souporcell to fluor based demuxing, pre-filtration
```{r}
#datasets I need to merge based on cell names. Here I want to not output all RNA reads
#but I want to output the Fluorescent reads, the soup assignments, the total RNA counts/cell, and
#the output of the multi calls on fluors.

#combine the dataframes by matching up cell barcodes
merge(soup_res_s1, multiDF_s1, by='row.names') -> mergedA_s1
#then have to make row names the row names again, manually
row.names(mergedA_s1)<-mergedA_s1$Row.names

#combine in next df
merge(mergedA_s1, FluorsDF_s1, by='row.names') -> mergedB_s1
#then have to make row names the row names again, manually
row.names(mergedB_s1)<-mergedB_s1$Row.names

#combine in next df
merge(mergedB_s1, RNAcountsDF_s1, by='row.names') -> mergedC_s1
#then, clean up DF columns from all extra row names cols
row.names(mergedC_s1)<-mergedC_s1$Row.names
#colnames(mergedC_s1)

#combine in next df
merge(mergedC_s1, FluorRNAcountsDF_s1, by='row.names') -> mergedD_s1
#then, clean up DF columns from all extra row names cols
row.names(mergedD_s1)<-mergedD_s1$Row.names
#colnames(mergedD_s1)

#change doublet calls (e.g. 0/1 to Doublet)
str_replace_all(mergedD_s1$soup_assign,"[:digit:]/[:digit:]", "Doublet") -> mergedD_s1$soup_assign

#now need to changes names in soup_assign to fluor proteins
table(mergedD_s1$multi_assign, mergedD_s1$soup_assign)

#so inspecting this it looks like mCherry expression is so low its likely worthless to
#worry too much about but given we know that this sample has: 
#7dpa 2 blastemas of 2 siblings(CAGGs:Venus)
#10dpa 3 blastemas  of 3 siblings (CAGGs:mCherry, B51)
#14dpa 3 blastemas  of 3 siblings (CAGGs:TFPnls, G48)
#there should be 2 Venus samples, and 3 Cherry and TFP in that case
#0 = Venus
#1 = TFP
#2 = Cherry (but likely not great due to low cherry reads)
#3 = Cherry (same as above)
#4 = ? due to low cell #, but based on sample #'s would guess TFP
#5 = Venus
#6 = TFP
#7 = Cherry

#based on table we will change to the above idents
str_replace_all(mergedD_s1$soup_assign,"[146]", "mTFP1") -> mergedD_s1$soup_assign
str_replace_all(mergedD_s1$soup_assign,"[05]", "Venus") -> mergedD_s1$soup_assign
str_replace_all(mergedD_s1$soup_assign,"[237]", "mCherryB51") -> mergedD_s1$soup_assign

#check to see if that workedout
table(mergedD_s1$multi_assign, mergedD_s1$soup_assign)

#merged$combined <-paste(merged$soup_assign, merged$cellplex_assign, sep = '_')
mergedD_s1$combined <-paste(mergedD_s1$soup_assign, mergedD_s1$multi_assign, sep = '_')
mergedD_s1$combined -> simple_s1
dim(mergedD_s1)
colnames(mergedD_s1)

dim(simple_s1)

table(simple_s1)

#remove extra 'Row.names columns' made by my inelegant combination of DFs by rowname
mergedD_s1[,c("soup_status","soup_assign","combined","multi_assign","RNACounts", "SumFluorCounts","Venus", "mTFP1", "mCherryB51")] -> df_XenowLabels

#this looks quite nice
ggplot(df_XenowLabels, aes(x=reorder(combined, combined, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))


```
## Now filter down to compare more granularly
```
```{r}
#saving of unfilterred dataframe (or just filterring done by CR)
write.table(df_XenowLabels, "/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/220803.Xenopus.fluor.soup.comparison.Unfiltered.txt", sep = '\t', col.names = T)

#filterring of Seurat object:
#setting a minimum and maximum reads counts. this is pretty typical except that we're raising the bar
#for minimum since we've seen this should help us just test accuracy on highest quality cells
#and then filterring for cells that have at least one read from a fluorophore marker gene
FeatureScatter(so.demux_s1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dim(so.demux_s1)
#13000 cells
so.demux_s1_filterred <- subset(so.demux_s1, subset = nCount_RNA > 5000 & nCount_RNA < 100000 & nCount_RNA > 0)
#down to 8530 cells
#takes me from around 13000 cells to 88530 or so.
FilterredCellBarcodes <- colnames(so.demux_s1_filterred)

#now want to filter our DF down to the same list of cells left from the above filterring, and save it.
Xeno.df.filterred <- df_XenowLabels[FilterredCellBarcodes,]

#now, save the filterredDF
write.table(Xeno.df.filterred,"/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/220803.Xeno.fluor.soup.comparison.Filterred.txt", sep = '\t', col.names = T)

```

##now compare filtered plots
First up, we want to plot summed fluorescent reads/cell versus avg 'usable calls' by multi on Fluorescent reads. binned every thousand cells, sorted by fluor reads/cell.
Next, we also may want to plot total reads/cell versus 'usable calls' by SOUP? and by Multi?
#in python for now

Then, we want Seurat tSNE plots of these filterred datasets, first colored by multi on Fluors.
Then colored by Soup assignments.
```{r}
#PCA: 
#normalizing:
so.demux_s1_filterred_Norm <- NormalizeData(so.demux_s1_filterred, normalization.method = "LogNormalize", scale.factor = 10000)

#find varialbe features
so.demux_s1_filterred_Norm <- FindVariableFeatures(so.demux_s1_filterred_Norm, selection.method = "vst", nfeatures = 2000)

#not sure if we need to first, normalize and scale all the data in the dataset?
#scale data. see seurat tutorial for why. says it's normal:
so.demux_s1_filterred_Norm <- ScaleData(so.demux_s1_filterred_Norm)
so.demux_s1_filterred_Norm <- RunPCA(so.demux_s1_filterred_Norm, features = VariableFeatures(object = so.demux_s1_filterred_Norm))

#is running umap needed?
so.demux_s1_filterred_Norm <- RunUMAP(so.demux_s1_filterred_Norm, dims = 1:25)

#coloring raw reads of each fluor
FeaturePlot(so.demux_s1_filterred_Norm, features = 'Venus') ->plot1
FeaturePlot(so.demux_s1_filterred_Norm, features = 'mCherryB51') ->plot2
FeaturePlot(so.demux_s1_filterred_Norm, features = 'mTFP1') ->plot3

#coloring by multi identity calls
DimPlot(so.demux_s1_filterred_Norm, group.by = 'MULTI_ID',cols = c("#E66100", "#D41159", "#40B0A6", "#AEB0C7",  "#FFC20A")) ->plot4

#coloring by soup calls
str_replace_all(so.demux_s1_filterred_Norm$assignment,"[:digit:]/[:digit:]", "Doublet") -> so.demux_s1_filterred_Norm$assignment
str_replace_all(so.demux_s1_filterred_Norm$assignment,"[146]", "mTFP1") -> so.demux_s1_filterred_Norm$assignment
str_replace_all(so.demux_s1_filterred_Norm$assignment,"[05]", "Venus") -> so.demux_s1_filterred_Norm$assignment
str_replace_all(so.demux_s1_filterred_Norm$assignment,"[237]", "mCherryB51") -> so.demux_s1_filterred_Norm$assignment

DimPlot(so.demux_s1_filterred_Norm, group.by = 'assignment', cols = c("#E66100", "#D41159", "#40B0A6", "#FFC20A")) ->plot5

#now plot
plot1 + plot2 + plot3 + plot4 + plot5
#save outputs when satisfied:


saveRDS(so.demux_s1_filterred_Norm, file = "/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/220803_XenoDemuxFilterredNormData.rds")
saveRDS(so.demux_s1, file = "/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/220803_XenoDemuxAllData.rds")
sessionInfo()

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_Venusplot.pdf") 
#or do same as above but png, takes up way less room
print(plot1)
dev.off()

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_mCherryplot.pdf") 
#or do same as above but png, takes up way less room
print(plot2)
dev.off()

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_TFPplot.pdf") 
#or do same as above but png, takes up way less room
print(plot3)
dev.off()

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_Multiplot.pdf") 
#or do same as above but png, takes up way less room
print(plot4)
dev.off()

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_Soupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot5)
dev.off()

#upset plot
so.demux_s1_filterred_Norm$combined <-paste(so.demux_s1_filterred_Norm$assignment, so.demux_s1_filterred_Norm$MULTI_ID, sep = '_')

ggplot(so.demux_s1_filterred_Norm@meta.data, aes(x=reorder(combined, combined, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90)) -> upset

pdf("/Users/nicholasleigh/Desktop/Lund/data/snp_paper/xenopus/early_timepoints/figures/220803_Xeno_Figure_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upset)
dev.off()

```

#now on to python code: https://colab.research.google.com/drive/1lO4ny8Uv9n1lPIbHmZKFPyxWI7gjmZPs?usp=sharing
