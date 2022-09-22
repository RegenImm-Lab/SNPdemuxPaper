
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

setwd('/Directory')
data_dir <- "/Directory/"

#Nick did a lot of filterring in R, here I'm just editing end of script to change the assignment names to allow direct comparisons and correlation analysis.

meta <- read.table('skip.remap.clusters.tsv', sep = '\t', header = T, row.names = 1)

table(meta$status)
#doublet    singlet unassigned 
#153       1786      19 

table(meta$assignment)
#0 0/1 0/2 0/3 0/4   1 1/0 1/2 1/3 1/4   2 2/0 2/1 2/3 2/4   3 3/0 3/1 3/2 3/4   4 4/0 4/1 4/2 
#273  22   9  14   9 370   5   7   6   6 368   4   7   1  17 255   3   7   8   3 521   5  17  13 
#4/3 
#12 


#the cell origin is sorted in the barcode info barcode-# with the # being the origin
#if it has a D it's a doublet and if it has S its a doublet but from the same animal

meta$barcodes <- rownames(meta)
separate(meta, barcodes, c("barcodes", "origin"), sep = "-",) -> meta

#going to clean things up a bit
meta -> cleaned.meta

#clean it up a bit by changing all #D in origin to doublet and all #/# in soup assignment to doublet
str_replace_all(cleaned.meta$assignment,"[:digit:]/[:digit:]", "Doublet") -> cleaned.meta$assignment
str_replace_all(cleaned.meta$origin,"[:digit:]D", "Doublet") -> cleaned.meta$origin

#remove lines with homotypic doublets since this is not a useful test for souporcell doublet detection
#troublet only detects doublets based on mixed genotypes
#cleaned.meta[cleaned.meta$origin != "[:digit:]S",] -> cleaned.metaTest
cleaned.meta[cleaned.meta$origin != "5S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "4S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "3S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "2S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "1S",] -> cleaned.meta

#match up highest matches
str_replace_all(cleaned.meta$assignment,"[0]", "5") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[2]", "Animal#Two") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[1]", "Animal#One") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[3]", "Animal#Three") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[4]", "Animal#Four") -> cleaned.meta$assignment

#now combine souporcell and ground truth to generate bar graph for an upset plot
cleaned.meta$combined_assignment <-paste(cleaned.meta$origin,cleaned.meta$assignment, sep = '_')

#upset plot
ggplot(cleaned.meta, aes(x=reorder(combined_assignment, combined_assignment, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->upsetPlot

pdf("Directory/2200917_GM_Figure_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upsetPlot)
dev.off()



#save for use in python
cleaned.meta[,c("status","assignment","combined_assignment", "origin")] -> df_Monkey

names(df_Monkey) <- c('soup_status', 'soup_assign', 'combined', 'origin')

write.table(df_Monkey, "/Directory/220917.GM.ground.truth.soup.comparison.txt", sep = '\t', col.names = T)
```
