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
Next, import data I need for spleen2 mapped to the dual species index
```{r}
#okay so this is gunna be quite a bit of loading
setwd('/Users/joecardiello/Documents/LeighLab/Data/Sp2_Dual_PleuroFlursNoto_ST_220718')

#load metadata from souporcell N=4 (dual species)
souporcell_calls_s2dual <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)
table(souporcell_calls_s2dual$assignment)

#will load in the 10X reads, and cellplex reads later so this isn't needed now.
#load in the 10X filterred read counts data
#Read10X(data.dir = "filtered_feature_bc_matrix/") -> rna_s2dual

#load in the cellplex data for each cell for the dual mapped data:
#load in data from cellranger making calls on cellplex
#cellplex <- read_csv('assignment_confidence_table.csv')
#as.data.frame(cellplex) -> df_cellplex
#rownames(df_cellplex) <- df_cellplex$Barcode
```

```{r}

#alternative option for souporcell read in is the souporcell n=3 run on spleen2 data
#mapped only to pleuro, so the 3 would be the two pleuro animals, and 1 for the combined noto
setwd('/Users/joecardiello/Documents/LeighLab/Data/Sp2_PleuroFluros_ST_220718')

#load metadata from souporcell N=3 (pleuor only)
souporcell_calls_s2pleuro <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)

table(souporcell_calls_s2pleuro$assignment)
#shows a lower doublet rate detected, and fewer total pleuro cells in each pleuro
#but also way fewer 'noto' type samples. so arguable which is better.
```


Quick simplification of souporcell info.
```{r}
#will simplify these a bit more
souporcell_calls_s2dual
souporcell_calls_s2dual[,c(1:2)] -> soup_res_s2dual

#so now need to see if soup_res when it calls #1 matching up most of the time with one of the CMO calls
dim(soup_res_s2dual)

#changing column names
names(soup_res_s2dual) <- c("soup_status", "soup_assign")
soup_res_s2dual
```

Loading in data matrixes from the raw outputs directory of the CR multi run. should
contain both the RNA reads from multi and the cellplex info.
Then filterring these down to cells in the CR counts, filterred, output.
```{r}
data_dir <- "/Users/joecardiello/Documents/LeighLab/Data/Sp2_Dual_PleuroFlursNoto_ST_220718/MultiOutputsRaw/raw_feature_bc_matrix"

sp2demux <- Read10X(data.dir = data_dir)

#rna counts 
sp2rna <- sp2demux[["Gene Expression"]]
#cellplex counts
sp2cellplex <- sp2demux[["Multiplexing Capture"]]

#just want the oligos we use. otherwise cellplex imports a list of around 10 CMOs that could potentially be used and makes a column for each. but we only used 3 labels in this experiment.
oligos <- c('CMO304', "CMO305", "CMO306")

#filter down cellplex
sp2cellplex[rownames(sp2cellplex) %in% oligos,] -> sp2small_cellplex

dim(sp2rna)
dim(sp2small_cellplex)

#try importing my filterred CR counts matrix as a way to filter down to cells that seem real:
setwd('/Users/joecardiello/Documents/LeighLab/Data/Sp2_Dual_PleuroFlursNoto_ST_220718')
Read10X(data.dir = "filtered_feature_bc_matrix/") -> Sp2counts
dim(Sp2counts)

#filter the rna matrix down to the list of cells from the filterred cells dataset
sp2joint.bcs <- intersect(colnames(sp2rna), colnames(Sp2counts))
length(sp2joint.bcs)

# Subset RNA and cellplex counts by joint cell barcodes, think this is already the case but good to do
sp2demux.umis <- sp2rna[, sp2joint.bcs]
sp2cellplex.counts <- as.matrix(sp2small_cellplex[, sp2joint.bcs])

# Confirm that the cellplex have the correct names and that filterring by real cells seems to have worked
rownames(sp2cellplex.counts)
dim(sp2demux.umis)
dim(sp2cellplex.counts)

#make seurat object with RNA reads and soup results as meta
#sp2so.demux <- CreateSeuratObject(counts = sp2demux.umis, project = "Round2_spl", assay = 'RNA')
sp2so.demux <- CreateSeuratObject(counts = sp2demux.umis, meta=soup_res_s2dual,project = "Round2_spl", assay = 'RNA')

#create a new assay for cellplex
sp2cellplex_assay <- CreateAssayObject(counts = sp2cellplex.counts)

# add this assay to the previously created Seurat object
sp2so.demux[["cellplex"]] <- sp2cellplex_assay

#this should output the read depth for each cell. and cellplex reads/cell
sp2so.demux$nCount_RNA -> sp2RNAcounts
sp2so.demux$nCount_cellplex -> sp2CellplexCounts
```

Now, try two different programs for calling cell identity on the cellplex info.
Try HTO and multi calling on cellplex. overall found multi works better for us.
```{r}
#normalize the lipid labelling data for these assays
sp2so.demux <- NormalizeData(sp2so.demux, assay = "cellplex", normalization.method = "CLR")

sp2so.demux <- HTODemux(sp2so.demux, assay = "cellplex", positive.quantile = 0.99)

# Global classification results
table(sp2so.demux$cellplex_classification)
table(sp2so.demux$cellplex_classification.global)
sp2so.demux$cellplex_classification-> sp2hto
#so.demux@meta.data$cellplex_classification

length(sp2so.demux$cellplex_classification.global)
length(sp2hto)

#good to try another method. Using Multi
sp2so.demux <- MULTIseqDemux(sp2so.demux, assay = 'cellplex', autoThresh = T)
#so.demux$MULTI_ID
length(sp2so.demux$MULTI_ID)
table(sp2so.demux$MULTI_ID)
#as.data.frame(so.demux$MULTI_ID) -> multiDF
sp2so.demux$MULTI_ID-> sp2multi
length(sp2multi)
#10474
```

Now just need to merge and save these output files for quick plotting of comparisons. will do this again after further calculations bellow.
After this, try working with these data sources to combine metadata into one df
to eventually enable comparing cell calls by cellplex vs soup etc
```{r}
#datasets to combine in one DF: souporcell calls, cellplex reads (full df), multi, and hto,
#plus the summed reads/cell for total RNA, and for cellplex labels
#to call these datasets, current names are: 
# soup_res_s2dual, sp2cellplex.counts, sp2multi,  #sp2hto, sp2RNAcounts, sp2CellplexCounts,

#make all releavnt datasets into DFs. reformatting or renaming columns as needed.
t(as.data.frame(sp2cellplex.counts)) -> sp2cellplex.RAWcountsDF
as.data.frame(sp2multi) ->sp2multiDF
names(sp2multiDF) <- c("multi_assign")
as.data.frame(sp2hto) ->sp2htoDF
names(sp2htoDF) <- c("hto_assign")
as.data.frame(sp2RNAcounts) ->sp2RNAcountsDF
names(sp2RNAcountsDF) <- c("RNACounts")
as.data.frame(sp2CellplexCounts) ->sp2CellplexCountsDF
names(sp2CellplexCountsDF) <- c("CellplexCounts")

#check out the  dfs/lists
dim(soup_res_s2dual)
dim(sp2cellplex.RAWcountsDF)
dim(sp2multiDF)
dim(sp2htoDF)
dim(sp2RNAcountsDF)
dim(sp2CellplexCountsDF)

#combine the dataframes by matching up cell barcodes
#repeatedly making sure I keep row labels with cell names. realize this is not an efficient way
merge(soup_res_s2dual, sp2cellplex.RAWcountsDF, by='row.names') -> mergedDFa
row.names(mergedDFa)<-mergedDFa$Row.names
mergedDFa$Barcode <-mergedDFa$Row.names

merge(mergedDFa, sp2multiDF, by='row.names') -> mergedDFb
row.names(mergedDFb)<-mergedDFb$Row.names

merge(mergedDFb, sp2htoDF, by='row.names') -> mergedDFc
row.names(mergedDFc)<-mergedDFc$Row.names

merge(mergedDFc, sp2RNAcountsDF, by='row.names') -> mergedDFd
row.names(mergedDFd)<-mergedDFd$Row.names

merge(mergedDFd, sp2CellplexCountsDF, by='row.names') -> mergedDFe
row.names(mergedDFe)<-mergedDFe$Row.names

#remove extra 'Row.names columns' made by my inelegant combination of DFs by rowname
mergedDFe[,c("Barcode","soup_status","soup_assign","hto_assign","multi_assign","RNACounts","CellplexCounts")] -> Sp2df_CellBarcodeLabels
dim(Sp2df_CellBarcodeLabels)


Sp2df_CellBarcodeLabels$combined <-paste(Sp2df_CellBarcodeLabels$soup_assign, Sp2df_CellBarcodeLabels$multi_assign, sep = '_')
Sp2df_CellBarcodeLabels$combined -> Sp2df_combinedLabels
dim(Sp2df_combinedLabels)
colnames(Sp2df_CellBarcodeLabels)

dim(Sp2df_combinedLabels)

table(Sp2df_combinedLabels)

#this looks quite nice
ggplot(Sp2df_CellBarcodeLabels, aes(x=reorder(combined, combined, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))


```

 To enable plotting/comparisons of calls, need to swap soup numeric calls to calls that match up with the cellplex names based on the above ggplot. Will do this again after further calculations bellow.
```{r}
#okay now that everything is in there...i think it makes most sense to change the soup_assign column
#to the most frequent assignemnt outside of the negative calls. which we can see in the above plot
#in relation to Multi, which is the cellplex assignment that seems to work best
#...that is 1=304, 0=305, 3=306, and 2=306 
#that way we should be able to make a plot that compares these when they all match
#make these assignments based on the above ggplot

#saving a copy of the merged DF for this name changing and table saving etc
Sp2dual.df <- Sp2df_CellBarcodeLabels
Sp2dual.df

#so flip out these in the soup_assign column
table(factor(Sp2dual.df$soup_assign))
table(factor(Sp2dual.df$multi_assign))
table(factor(Sp2dual.df$hto_assign))
#see order to be replaced

Sp2dual.df$NEWsoup_assign <- factor(Sp2dual.df$soup_assign, labels = c('CMO305','Doublet', 'Doublet', 'Doublet','CMO304', 'Doublet', 'Doublet', 'Doublet','CMO306', 'Doublet', 'Doublet', 'Doublet','CMO306', 'Doublet', 'Doublet', 'Doublet'))
table(factor(Sp2dual.df$NEWsoup_assign))

#need to remove "unassingned" and blanks from cellplex  to call Negative. and repalce multiplet with doublet
str_replace_all(Sp2dual.df$soup_status, 'unassigned', 'Negative') -> test_sp2dual
str_replace_all(test_sp2dual, 'doublet', 'Doublet') -> test_sp2dual
#str_replace_all(test2, 'Blanks', 'Negative') -> test3
#table(test3)

Sp2dual.df$soup_status <- test_sp2dual
table(Sp2dual.df$NEWsoup_assign)
```

Next, due to the constraints in the read depth/fluor read depth vs accuracy/cell calls samples,
we want to make a copy of the DF and Seurat object, 
then filter the DF, and the Seurat object by a lowest total read depth
Want to save a file of the Dataframe before and after Filtration, to enable plotting of both
```{r}
#saving of unfilterred dataframe (or just filterring done by CR)
write.table(Sp2dual.df, "/Users/joecardiello/Documents/LeighLab/Data/220804.Figure4.spleen2.cellplex.soup.comparison.Unfilterred.txt", sep = '\t', col.names = T)

#filterring of Seurat object:
#setting a minimum and maximum reads counts. this is pretty typical except that we're raising the bar for minimum since we've seen this should help us just test accuracy on highest quality cells and then filterring for cells that have at least one read from a fluorophore marker gene
sp2so.demuxFilterred <- subset(sp2so.demux, subset = nCount_RNA > 5000 & nCount_RNA < 40000)

#takes me from around 10474 to 4043 cells or so.
sp2so.demuxFilterredBarcodes <- colnames(sp2so.demuxFilterred)

#now want to filter our DF down to the same list of cells left from the above filterring, and save it.
Sp2dual.df.filterred <- Sp2dual.df[sp2so.demuxFilterredBarcodes,]

#now, save the filterredDF
write.table(Sp2dual.df.filterred, "/Users/joecardiello/Documents/LeighLab/Data/220804.Figure4.spleen2.cellplex.soup.comparison.Filterred.txt", sep = '\t', col.names = T)
```

First up, we want to plot summed fluorescent reads/cell versus avg 'usable calls' by multi on Fluorescent reads. binned every thousand cells, sorted by fluor reads/cell.
Next, we also may want to plot total reads/cell versus 'usable calls' by SOUP? and by Multi?
#in python for now

Then, we want Seurat tSNE plots of these filterred datasets, first colored by multi on Fluors.
Then colored by Soup assignments.
```{r}
#PCA: 
#normalizing:
sp2so.demuxFilterred_norm <- NormalizeData(sp2so.demuxFilterred, normalization.method = "LogNormalize", scale.factor = 10000)

#Kept getting errors with memory on the scaling. found suggestions to just scale the
#variable features, to reduce memory. #this does take forever. 5-10 min?
sp2so.demuxFilterred_norm <- FindVariableFeatures(sp2so.demuxFilterred_norm, selection.method = "vst", nfeatures = 2000)

#not sure if we need to first, normalize and scale all the data in the dataset?
#scale data. see seurat tutorial for why. says it's normal:
sp2so.demuxFilterred_norm <- ScaleData(sp2so.demuxFilterred_norm)
sp2so.demuxFilterred_norm <- RunPCA(sp2so.demuxFilterred_norm, features = VariableFeatures(object = sp2so.demuxFilterred_norm))

#is running umap needed?
sp2so.demuxFilterred_norm <- RunUMAP(sp2so.demuxFilterred_norm, dims = 1:25)
```
Plotting and saving plots:
```{r}
#coloring raw reads of some genes of interest
FeaturePlot(sp2so.demuxFilterred_norm, features = 'CMO304') ->plot4_1
FeaturePlot(sp2so.demuxFilterred_norm, features = 'CMO305')->plot4_2
FeaturePlot(sp2so.demuxFilterred_norm, features = 'CMO306')->plot4_3

#coloring by cellplex multi calls
DimPlot(sp2so.demuxFilterred_norm, group.by = 'MULTI_ID') ->plot4_4

#coloring by soup calls
DimPlot(sp2so.demuxFilterred_norm, group.by = 'soup_assign') ->plot4_5

#save outputs when satisfied:
saveRDS(sp2so.demuxFilterred_norm, file = "/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Spleen2CellplexDualDemuxFilterredNormData.rds")
#saveRDS(so.demux_s2_dual, file = "/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Spleen2DualDemuxAllData.rds")
sessionInfo()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_Cellplex1plot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_1)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_Cellplex2plot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_2)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_Cellplex3plot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_3)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_MultiCellplexColorsUmapplot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_4)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_SoupUmapPlot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_5)
dev.off()

table(factor(sp2so.demuxFilterred_norm$MULTI_ID))
table(factor(sp2so.demuxFilterred_norm$soup_assign))
#see order to be replaced
#...that is 1=304, 0=305, 3=306, and 2=306 
sp2so.demuxFilterred_norm$soup_assign <- factor(sp2so.demuxFilterred_norm$soup_assign, labels = c('CMO305','Doublet', 'Doublet', 'Doublet','CMO304', 'Doublet', 'Doublet', 'Doublet','CMO306', 'Doublet', 'Doublet','CMO306', 'Doublet', 'Doublet', 'Doublet'))
table(factor(sp2so.demuxFilterred_norm$soup_assign))
#coloring by soup calls
DimPlot(sp2so.demuxFilterred_norm, group.by = 'soup_assign') -> plot4_6
#DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'assignment',cols = c("blue", "purple", "red", "green")) -> plot6
pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_SimplifiedSoupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot4_6)
dev.off()

```


Next, similar to the ggplot near the start, we want an upset plots of cell assignments by these two methods.
```{r}
Sp2dual.df.filterred$combined2 <-paste(Sp2dual.df.filterred$multi_assign, Sp2dual.df.filterred$NEWsoup_assign, sep = '_')
ggplot(Sp2dual.df.filterred, aes(x=reorder(combined2, combined2, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->plot4_7

#table(simple_s1)

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig4_Sp2_Cellplex/220804_Figure4_Upset.pdf") 
#or do same as above but png, takes up way less room
print(plot4_7)
dev.off()

```
