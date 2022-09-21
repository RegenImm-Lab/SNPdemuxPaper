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

Read10X(data.dir = "filtered_feature_bc_matrix/") -> rna_s2dual

#load metadata from souporcell N=4 (dual species)
souporcell_calls_s2dual <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)

table(souporcell_calls_s2dual$assignment)

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

Making a Seurat object and DF of RNA total counts/cell
```{r}
#make seurat object with the gene expression counts,
#and with souporcell results as the metadata
so.demux_s2_dual <- CreateSeuratObject(counts = rna_s2dual, meta=soup_res_s2dual, project = "Round2_spl", assay = 'RNA')
#so.demux_s2_dual <- CreateSeuratObject(counts = rna_s2dual, meta=souporcell_calls_s2dual, project = "Round2_spl", assay = 'RNA')

#this should output the read depth for each cell. but I need to confirm that this automated
#function is working correctly
so.demux_s2_dual$nCount_RNA -> RNAcounts_s2_dual
as.data.frame(RNAcounts_s2_dual) -> RNAcounts_s2_dual_DF
#changing column names in the DF
names(RNAcounts_s2_dual_DF) <- c("RNACounts")

```
Here I'll try out some methods to analyze this CR dataset
Trying out part of Nick's method for barnyard analysis of these data:
```{r}
#use seurat subset to filter only cells with more than 5000 counts
#I will do this later downstream.
#filtered <- subset(r2_Seu, subset = nCount_RNA > 5000)
dim(so.demux_s2_dual)
#[1] 1123219   10474


#make a matrix with just counts info
GEX <- so.demux_s2_dual@assays$RNA@counts

#now see if we can look at barnyard plot
#using code from https://www.kallistobus.tools/tutorials/kb_species_mixing/r/kb_mixed_species_10x_v2/
gene_species <- ifelse(str_detect(rownames(GEX), "^NotoST"), "Noto", "Pleuro")
Pleuro_inds <- gene_species == "Pleuro"
Noto_inds <- gene_species == "Noto"
# mark cells as Pleuro or Noto
cell_species <- tibble(n_Pleuro_umi = Matrix::colSums(GEX[Pleuro_inds,]),
                       n_Noto_umi = Matrix::colSums(GEX[Noto_inds,]),
                       tot_umi = Matrix::colSums(GEX),
                       prop_Pleuro = n_Pleuro_umi / tot_umi,
                       prop_Noto = n_Noto_umi / tot_umi, 
                       barcodes = colnames(GEX))

# Classify species based on proportion of UMI, with cutoff of 90% or adjusted to what you want
cell_species <- cell_species %>% 
  mutate(species = case_when(
    prop_Pleuro> 0.9 ~ "Pleuro",
    prop_Noto > 0.6 ~ "Noto",
    TRUE ~ "Doublet"
  ))

#barnyard plot
ggplot(cell_species, aes(n_Noto_umi, n_Pleuro_umi, color = species)) +
  geom_point()

#print my barnyard plot after filterring:
ggplot(cell_species, aes(n_Noto_umi, n_Pleuro_umi, color = species)) +
  geom_point() -> BarnPlotSp2Dual_unfilterred

pdf("/Users/joecardiello/Documents/LeighLab/Data/220809_DemuxFiguresV2/Fig2_Sp2_Barnyard/220907_Figure2_BarnyardPlotUnfilterred.pdf") 
#or do same as above but png, takes up way less room
print(BarnPlotSp2Dual_unfilterred)
dev.off()

#tibble to dataframe for loading as metadata into seurat object
as.data.frame(cell_species) -> df_cell_sp

#make colnames barcodes for metadata addition
rownames(df_cell_sp) <- df_cell_sp$barcodes

#to Seurat object that already has RNA counts and souporcell results as metadata
#now add the barnyard assignments as second set of metadata:
so.demux_s2_dual<-AddMetaData(so.demux_s2_dual,metadata = df_cell_sp )

#table(simple_s1)



```


Merge first datasets, and use ggplot to determine which souporcell label likely matches with
each fluorescent label.
```{r}
#datasets I need to merge based on cell names. Here I want to not output all RNA reads
#but I want to output the barnyard assignments/counts, the soup assignments, and the total RNA counts/cell
#combine: soup_res_s2dual,  df_cell_sp,RNAcounts_s2_dual_DF

#combine the dataframes by matching up cell barcodes
merge(soup_res_s2dual, df_cell_sp, by='row.names') -> mergedA_s2_dual
#then have to make row names the row names again, manually
row.names(mergedA_s2_dual)<-mergedA_s2_dual$Row.names

#combine in next df
merge(mergedA_s2_dual, RNAcounts_s2_dual_DF, by='row.names') -> mergedA_s2_dual2
#then have to make row names the row names again, manually
row.names(mergedA_s2_dual2)<-mergedA_s2_dual2$Row.names


mergedA_s2_dual2$combined <-paste(mergedA_s2_dual2$soup_assign, mergedA_s2_dual2$species, sep = '_')
mergedA_s2_dual2$combined -> simple_s2dual
dim(simple_s2dual)
colnames(mergedA_s2_dual2)

dim(simple_s2dual)

table(simple_s2dual)

#remove extra 'Row.names columns' made by my inelegant combination of DFs by rowname
mergedA_s2_dual2[,c("soup_status","soup_assign","combined","species","RNACounts", "tot_umi","n_Pleuro_umi", "n_Noto_umi", "prop_Pleuro","prop_Noto")] -> df_Spleen2DualwLabels

#this looks quite nice
ggplot(df_Spleen2DualwLabels, aes(x=reorder(combined, combined, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))

```

 To enable plotting/comparisons of calls, need to swap soup numeric calls to calls that match up with the cellplex names based on the above ggplot. Will do this again after further calculations bellow.
```{r}
#okay now that everything is in there...i think it makes most sense to change the soup_assign column
#to the most frequent assignemnt outside of the negative calls
#...that is 1=pleuro, 0=pleuro, 3=noto, and 2=noto (or mixed)
#that way we should be able to make a plot that compares these when they all match
#make these assignments based on the above ggplot

#saving a copy of the merged DF for this name changing and table saving etc
Spleen2dual.df <- df_Spleen2DualwLabels
Spleen2dual.df

#so flip out these in the soup_assign column
table(factor(Spleen2dual.df$soup_assign))
table(factor(Spleen2dual.df$species))
#see order to be replaced

Spleen2dual.df$NEWsoup_assign <- factor(Spleen2dual.df$soup_assign, labels = c('Pleuro','Doublet', 'Doublet', 'Doublet','Pleuro', 'Doublet', 'Doublet', 'Doublet','Noto', 'Doublet', 'Doublet', 'Doublet','Noto', 'Doublet', 'Doublet', 'Doublet'))
table(factor(Spleen2dual.df$NEWsoup_assign))
#cleaned.df

#need to remove "unassingned" and blanks from cellplex adn repalce multiplet with doublet

str_replace_all(Spleen2dual.df$soup_status, 'unassigned', 'Negative') -> test_s2dual
str_replace_all(test_s2dual, 'doublet', 'Doublet') -> test2_s2dual
#str_replace_all(test2, 'Blanks', 'Negative') -> test3
#table(test3)

Spleen2dual.df$soup_status <- test2_s2dual

table(Spleen2dual.df$NEWsoup_assign)


#save DF if I just want to compare cellplex and souporcell
```

Next, due to the constraints in the read depth/fluor read depth vs accuracy/cell calls samples,
we want to make a copy of the DF and Seurat object, 
then filter the DF, and the Seurat object for cells with 0 Fluorescent Reads.
And by a lowest total read depth??
Want to save a file of the Dataframe after Filtration, to enable plotting of both
```{r}
#saving of unfilterred dataframe (or just filterring done by CR)
write.table(Spleen2dual.df, "/Users/joecardiello/Documents/LeighLab/Data/220729.Figure2.spleen2.barnyard.soup.comparison.Unfilterred.txt", sep = '\t', col.names = T)

#filterring of Seurat object:
#setting a minimum and maximum reads counts. this is pretty typical except that we're raising the bar
#for minimum since we've seen this should help us just test accuracy on highest quality cells
#and then filterring for cells that have at least one read from a fluorophore marker gene
so.demux_s2_dual_filterred <- subset(so.demux_s2_dual, subset = nCount_RNA > 5000 & nCount_RNA < 40000)
#losing my 6 columns of metadata included in the barnyard experiment here for some reason.
#think this isn't actually running...

#takes me from around 10474 to 4043 cells or so.
FilterredCellBarcodesSp2Dual <- colnames(so.demux_s2_dual_filterred)

#now want to filter our DF down to the same list of cells left from the above filterring, and save it.
Spleen2dual.df.filterred <- Spleen2dual.df[FilterredCellBarcodesSp2Dual,]

#now, save the filterredDF
write.table(Spleen2dual.df.filterred, "/Users/joecardiello/Documents/LeighLab/Data/220729.Figure2.spleen2.barnyard.soup.comparison.Filterred.txt", sep = '\t', col.names = T)
table(Spleen2dual.df.filterred$species)

```

We want Seurat tSNE plots of these filterred datasets, first colored by barnyard
Then colored by Soup assignments.
```{r}
#PCA: 
#normalizing:
so.demux_s2_dual_filterred_norm <- NormalizeData(so.demux_s2_dual_filterred, normalization.method = "LogNormalize", scale.factor = 10000)

#Kept getting errors with memory on the scaling. found suggestions to just scale the
#variable features, to reduce memory. #this does take forever. 5-10 min?
so.demux_s2_dual_filterred_norm <- FindVariableFeatures(so.demux_s2_dual_filterred_norm, selection.method = "vst", nfeatures = 2000)

#not sure if we need to first, normalize and scale all the data in the dataset?
#scale data. see seurat tutorial for why. says it's normal:
so.demux_s2_dual_filterred_norm <- ScaleData(so.demux_s2_dual_filterred_norm)
so.demux_s2_dual_filterred_norm <- RunPCA(so.demux_s2_dual_filterred_norm, features = VariableFeatures(object = so.demux_s2_dual_filterred_norm))

#is running umap needed?
so.demux_s2_dual_filterred_norm <- RunUMAP(so.demux_s2_dual_filterred_norm, dims = 1:25)
```
Plotting and saving plots:
```{r}
#coloring raw reads of some genes of interest
FeaturePlot(so.demux_s2_dual_filterred_norm, features = 'PleuroST-GFP') ->plot1

FeaturePlot(so.demux_s2_dual_filterred_norm, features = 'NotoST---43158')->plot2
FeaturePlot(so.demux_s2_dual_filterred_norm, features = 'NotoST---30451')->plot3

#FeaturePlot(so.demux_s1_filterred_Norm, features = 'CHERRY') ->plot2
#FeaturePlot(so.demux_s1_filterred_Norm, features = 'EBFP') ->plot3

#coloring by barnyard species calls
DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'species') ->plot4

#coloring by soup calls
DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'soup_assign')
DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'soup_assign') ->plot5

#save outputs when satisfied:

saveRDS(so.demux_s2_dual_filterred_norm, file = "/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Spleen2DualDemuxFilterredNormData.rds")
saveRDS(so.demux_s2_dual, file = "/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Spleen2DualDemuxAllData.rds")
sessionInfo()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_GFPplot.pdf") 
#or do same as above but png, takes up way less room
print(plot1)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_NotoGene1yplot.pdf") 
#or do same as above but png, takes up way less room
print(plot2)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_NotoGene2plot.pdf") 
#or do same as above but png, takes up way less room
print(plot3)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_BanyardColorsUmapplot.pdf") 
#or do same as above but png, takes up way less room
print(plot4)
dev.off()

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_Soupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot5)
dev.off()

table(factor(so.demux_s2_dual_filterred_norm$species))
table(factor(so.demux_s2_dual_filterred_norm$soup_assign))
#see order to be replaced
#...that is 1=pleuro, 0=pleuro, 3=noto, and 2=noto (or mixed)
#labels = c('pleuro','Doublet', 'Doublet', 'Doublet','pleuro', 'Doublet', 'Doublet', 'Doublet','noto', 'Doublet', 'Doublet', 'Doublet','noto', 'Doublet', 'Doublet', 'Doublet')
so.demux_s2_dual_filterred_norm$soup_assign <- factor(so.demux_s2_dual_filterred_norm$soup_assign, labels = c('Pleuro','Doublet', 'Doublet', 'Doublet','Pleuro', 'Doublet', 'Doublet', 'Doublet','Noto', 'Doublet', 'Doublet','Noto', 'Doublet', 'Doublet', 'Doublet'))
table(factor(so.demux_s2_dual_filterred_norm$soup_assign))
#coloring by soup calls
DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'soup_assign') -> plot6
#DimPlot(so.demux_s2_dual_filterred_norm, group.by = 'assignment',cols = c("blue", "purple", "red", "green")) -> plot6
pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_SimplifiedSoupplot.pdf") 
#or do same as above but png, takes up way less room
print(plot6)
dev.off()

```


Next, similar to the ggplot near the start, we want an upset plots of cell assignments by these two methods.
```{r}
Spleen2dual.df.filterred$combined2 <-paste(Spleen2dual.df.filterred$species, Spleen2dual.df.filterred$NEWsoup_assign, sep = '_')
ggplot(Spleen2dual.df.filterred, aes(x=reorder(combined2, combined2, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->upsetPlotSp2Dual

#table(simple_s1)

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220729_Figure2_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upsetPlotSp2Dual)
dev.off()

#print my barnyard plot after filterring:

ggplot(Spleen2dual.df.filterred, aes(n_Noto_umi, n_Pleuro_umi, color = species)) +
  geom_point() ->BarnPlotSp2Dual

pdf("/Users/joecardiello/Documents/LeighLab/Data/220802_DemuxFigures/Fig2_Sp2_Barnyard/220801_Figure2_BarnyardPlot.pdf") 
#or do same as above but png, takes up way less room
print(BarnPlotSp2Dual)
dev.off()
```
