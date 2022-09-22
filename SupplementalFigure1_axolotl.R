Code from Nick for analyzing synthetically pooled axolotl data
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
data_dir <- "/Directory"

#meta <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)
meta <- read.table('axo.clusters.tsv', sep = '\t', header = T)

readnumb <- read.table('axo.snp_containing_reads_per_barcode', sep = '',col.names=c("SNPreads","barcode"))
separate(readnumb, barcode, sep = "-", c("barcodes", "originAgain")) -> readnumb

table(meta$status)
#doublet    singlet unassigned 
#334      10085        342 

table(meta$assignment)
#0  0/1  0/2    1  1/0  1/2    2  2/0  2/1 
#3286   85  126 1875   66   46 5068   94  115 

#from web summaries cell numbers for 
#D_1 = 2,193
#L_1 = 3,628
#M_1 = 6,016

#will keep unassinged 
#the cell origin is sorted in the barcode info barcode-# with the # being the origin
#if it has a D it's a doublet and if it has S its a doublet but from the same animal
separate(meta, barcode, sep = "-", c("barcodes", "origin")) -> meta
meta <- merge(meta,readnumb, by = "barcodes")

#going to clean things up a bit
meta -> cleaned.meta

#clean it up a bit by changing all #D in origin to doublet and all #/# in soup assignment to doublet
str_replace_all(cleaned.meta$assignment,"[:digit:]/[:digit:]", "Doublet") -> cleaned.meta$assignment
str_replace_all(cleaned.meta$origin,"[:digit:]D", "Doublet") -> cleaned.meta$origin

#remove lines with homotypic doublets since this is not a useful test for souporcell doublet detection
#troublet only detects doublets based on different genotypes
#cleaned.meta[cleaned.meta$origin != "[:digit:]S",] -> cleaned.metaTest
cleaned.meta[cleaned.meta$origin != "3S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "2S",] -> cleaned.meta
cleaned.meta[cleaned.meta$origin != "1S",] -> cleaned.meta

#here we're replacing the souporcell assignments with values matching the animal origin assigments
#that were found to correlate with them via an upset plot.
#str_replace_all(cleaned.meta$assignment,"[2]", "Animal#3") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[0]", "Animal#2") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[1]", "Animal#1") -> cleaned.meta$assignment
str_replace_all(cleaned.meta$assignment,"[2]", "3") -> cleaned.meta$assignment
str_replace_all(cleaned.meta$assignment,"[0]", "2") -> cleaned.meta$assignment
str_replace_all(cleaned.meta$assignment,"[1]", "1") -> cleaned.meta$assignment


#now combine souporcell and ground truth to generate bar graph for an upset plot
cleaned.meta$combined_assignment <-paste(cleaned.meta$origin,cleaned.meta$assignment, sep = '_')

#upset plot
ggplot(cleaned.meta, aes(x=reorder(combined_assignment, combined_assignment, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->upsetPlot

pdf("/Directory/220825_Axo_Figure_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upsetPlot)
dev.off()

#save for use in python
cleaned.meta[,c("status","assignment","combined_assignment", "origin","SNPreads")] -> df_Axo

names(df_Axo) <- c('soup_status', 'soup_assign', 'combined', 'origin',"SNPreads")

write.table(df_Axo, "/Directory/220825.Axo.ground.truth.soup.comparison.txt", sep = '\t', col.names = T)

```
