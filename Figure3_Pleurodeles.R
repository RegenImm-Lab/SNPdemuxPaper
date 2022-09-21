
my code adaptations from Nick's bellow. First, get programs we need
```{r}
library(Seurat)
library(tximport)
library(dplyr)
library(stringr)
library(ggplot2)
library(BUSpaRse)
library(Seurat)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(tidyr)
```
Next, import data I need for spleen1
```{r}

setwd('/Users/joecardiello/Documents/LeighLab/Data/Sp1_PleuroFluros_ST_220718')

Read10X(data.dir = "filtered_feature_bc_matrix/") -> rna_s1

#load metadata from souporcell
souporcell_calls_s1 <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)


```


Quick simplification of souporcell info.
```{r}
#will simplify these a bit more
souporcell_calls_s1
souporcell_calls_s1[,c(1:2)] -> soup_res_s1

#so now need to see if soup_res when it calls #1 matching up most of the time with one of the CMO calls
dim(soup_res_s1)

#changing column names
names(soup_res_s1) <- c("soup_status", "soup_assign")

```

Here I'll try out some methods to analyze this CR dataset 
Looking at the fluorophore counts as a metric for IDing cell identity.
```{r}
#try to print the subset of this gene expression data that's the fluors:
#rna['EBFP',]
oligos <- c('EBFP', "GFP", "CHERRY")

#make DF of Fluor counts similar to Cellplex strategy
rna_s1[rownames(rna_s1) %in% oligos,] -> FluorsGenes
as.data.frame(FluorsGenes) -> FluorsDFa
t(FluorsDFa) -> FluorsDF_s1

#make seurat object with the gene expression counts,
#and with souporcell results as the metadata
#so.demux_s1 <- CreateSeuratObject(counts = rna_s1, project = "Round1_spl", assay = 'RNA')
so.demux_s1 <- CreateSeuratObject(counts = rna_s1, meta=souporcell_calls_s1, project = "Round1_spl", assay = 'RNA')

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



```

Now, try two different programs for calling cell identity on the cellplex info.
Try applying this to the calls on fluorophore counts.
```{r}
#normalize the fluorescent labelling data for these assays similarly to cellplex
so.demux_s1 <- NormalizeData(so.demux_s1, assay = "fluors", normalization.method = "CLR")

#try HTO
#HTO is failing for fluorophores
#so.demux_s1 <- HTODemux(so.demux_s1, assay = "fluors", positive.quantile = 0.99)
# Global classification results
#table(so.demux$fluors_classification)
#table(so.demux$fluors_classification.global)
#so.demux$fluors_classification-> hto
#so.demux@meta.data$cellplex_classification

#length(so.demux$cellplex_classification.global)
#length(hto)

#we've found multi works better for us than HTO
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


Merge first datasets, and use ggplot to determine which souporcell label likely matches with
each fluorescent label.
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

#merged$combined <-paste(merged$soup_assign, merged$cellplex_assign, sep = '_')
mergedD_s1$combined <-paste(mergedD_s1$soup_assign, mergedD_s1$multi_assign, sep = '_')
mergedD_s1$combined -> simple_s1
dim(mergedD_s1)
colnames(mergedD_s1)

dim(simple_s1)

table(simple_s1)

#remove extra 'Row.names columns' made by my inelegant combination of DFs by rowname
mergedD_s1[,c("soup_status","soup_assign","combined","multi_assign","RNACounts", "SumFluorCounts","GFP", "CHERRY", "EBFP")] -> df_Spleen1wLabels

#this looks quite nice
ggplot(df_Spleen1wLabels, aes(x=reorder(combined, combined, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))


```

 To enable plotting/comparisons of calls, need to swap soup numeric calls to calls that match up with the cellplex names based on the above ggplot. Will do this again after further calculations bellow.
```{r}
#okay now that everything is in there...i think it makes most sense to change the soup_assign column
#to the most frequent assignemnt outside of the negative calls
#...that is 2 = GFP, 1=CHERRY and 0=EBFP
#that way we should be able to make a plot that compares these when they all match
#make these assignments based on the above ggplot

#saving a copy of the merged DF for this name changing and table saving etc
Spleen1.df <- df_Spleen1wLabels
Spleen1.df

#so flip out these in the soup_assign column
table(factor(Spleen1.df$soup_assign))
table(factor(Spleen1.df$multi_assign))
#see order to be replaced

Spleen1.df$NEWsoup_assign <- factor(Spleen1.df$soup_assign, labels = c('EBFP','Doublet', 'Doublet', 'CHERRY', 'Doublet', 'Doublet','GFP', 'Doublet', 'Doublet'))
table(factor(Spleen1.df$NEWsoup_assign))
#cleaned.df

#need to remove "unassingned" and blanks from cellplex adn repalce multiplet with doublet

str_replace_all(Spleen1.df$soup_status, 'unassigned', 'Negative') -> test_s1
str_replace_all(test_s1, 'doublet', 'Doublet') -> test2_s1
#str_replace_all(test2, 'Blanks', 'Negative') -> test3
#table(test3)

Spleen1.df$soup_status <- test2_s1

table(Spleen1.df$NEWsoup_assign)


#save DF if I just want to compare cellplex and souporcell
```

Next, due to the constraints in the read depth/fluor read depth vs accuracy/cell calls samples,
we want to make a copy of the DF and Seurat object, 
then filter the DF, and the Seurat object for cells with 0 Fluorescent Reads.
And by a lowest total read depth??
Want to save a file of the Dataframe after Filtration, to enable plotting of both
```{r}
#saving of unfilterred dataframe (or just filterring done by CR)
write.table(Spleen1.df, "/Users/joecardiello/Documents/LeighLab/Data/220721.spleen1.fluor.soup.comparison.Unfilterred.txt", sep = '\t', col.names = T)

#filterring of Seurat object:
#setting a minimum and maximum reads counts. this is pretty typical except that we're raising the bar
#for minimum since we've seen this should help us just test accuracy on highest quality cells
#and then filterring for cells that have at least one read from a fluorophore marker gene
so.demux_s1_filterred <- subset(so.demux_s1, subset = nCount_RNA > 5000 & nCount_RNA < 40000 & nCount_RNA > 0)
#takes me from around 25,700 cells to 5,800 or so.
FilterredCellBarcodes <- colnames(so.demux_s1_filterred)

#now want to filter our DF down to the same list of cells left from the above filterring, and save it.
Spleen1.df.filterred <- Spleen1.df[FilterredCellBarcodes,]

#now, save the filterredDF
write.table(Spleen1.df.filterred, "/Users/joecardiello/Documents/LeighLab/Data/220721.spleen1.fluor.soup.comparison.Filterred.txt", sep = '\t', col.names = T)

```


First up, we want to plot summed fluorescent reads/cell versus avg 'usable calls' by multi on Fluorescent reads. binned every thousand cells, sorted by fluor reads/cell.
Next, we also may want to plot total reads/cell versus 'usable calls' by SOUP? and by Multi?
#in python for now

Then, we want Seurat tSNE plots of these filterred datasets, first colored by multi on Fluors.
Then colored by Soup assignments.
```{r}
#PCA: 
#normalizing:
so.demux_s1_filterred_Norm <- NormalizeData(so.demux_s1_filterred, normalization.method = "LogNormalize", scale.factor = 10000)

#Kept getting errors with memory on the scaling. found suggestions to just scale the
#variable features, to reduce memory. #this does take forever. 5-10 min?
so.demux_s1_filterred_Norm <- FindVariableFeatures(so.demux_s1_filterred_Norm, selection.method = "vst", nfeatures = 2000)

#not sure if we need to first, normalize and scale all the data in the dataset?
#scale data. see seurat tutorial for why. says it's normal:
so.demux_s1_filterred_Norm <- ScaleData(so.demux_s1_filterred_Norm)
so.demux_s1_filterred_Norm <- RunPCA(so.demux_s1_filterred_Norm, features = VariableFeatures(object = so.demux_s1_filterred_Norm))

#is running umap needed?
so.demux_s1_filterred_Norm <- RunUMAP(so.demux_s1_filterred_Norm, dims = 1:25)

#coloring raw reads of each fluor
FeaturePlot(so.demux_s1_filterred_Norm, features = 'GFP') ->plot1
FeaturePlot(so.demux_s1_filterred_Norm, features = 'CHERRY') ->plot2
FeaturePlot(so.demux_s1_filterred_Norm, features = 'EBFP') ->plot3

#coloring by multi identity calls
DimPlot(so.demux_s1_filterred_Norm, group.by = 'MULTI_ID',cols = c("red", "purple", "blue", "green","pink")) ->plot4

#coloring by soup calls
DimPlot(so.demux_s1_filterred_Norm, group.by = 'assignment') ->plot5

#save outputs when satisfied:

saveRDS(so.demux_s1_filterred_Norm, file = "/Users/joecardiello/Documents/LeighLab/Data/220722_DemuxRobjects/220722_Spleen1DemuxFilterredNormData.rds")
saveRDS(so.demux_s1, file = "/Users/joecardiello/Documents/LeighLab/Data/220722_DemuxRobjects/220722_Spleen1DemuxAllData.rds")
sessionInfo()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_GFPplot.pdf") 
#or do same as above but png, takes up way less room
print(plot1)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_Cherryplot.pdf") 
#or do same as above but png, takes up way less room
print(plot2)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_EBFPplot.pdf") 
#or do same as above but png, takes up way less room
print(plot3)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors//220722_Figure3_Multiplot.pdf") 
#or do same as above but png, takes up way less room
print(plot4)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_Soupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot5)
dev.off()

table(factor(so.demux_s1_filterred_Norm$MULTI_ID))
table(factor(so.demux_s1_filterred_Norm$assignment))
#see order to be replaced
#2 = GFP, 1=CHERRY and 0=EBFP
so.demux_s1_filterred_Norm$assignment <- factor(so.demux_s1_filterred_Norm$assignment, labels = c('EBFP','Doublet', 'Doublet', 'CHERRY', 'Doublet', 'Doublet','GFP', 'Doublet', 'Doublet'))
table(factor(so.demux_s1_filterred_Norm$assignment))
#coloring by soup calls
DimPlot(so.demux_s1_filterred_Norm, group.by = 'assignment',cols = c("blue", "purple", "red", "green")) -> plot6
pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_SimplifiedSoupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot6)
dev.off()

```


Next, similar to the ggplot near the start, we want an upset plot of cell assignments by these two methods.
```{r}
Spleen1.df.filterred$combined2 <-paste(Spleen1.df.filterred$multi_assign, Spleen1.df.filterred$NEWsoup_assign, sep = '_')
ggplot(Spleen1.df.filterred, aes(x=reorder(combined2, combined2, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->upsetPlot


table(simple_s1)

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig3_Sp1_Fluors/220722_Figure3_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upsetPlot)
dev.off()
```
