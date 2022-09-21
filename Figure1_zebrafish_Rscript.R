
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

setwd('/Users/joecardiello/Documents/LeighLab/Data/220827_ZFdata')
data_dir <- "/Users/joecardiello/Documents/LeighLab/Data/"

#meta <- read.table('clusters.tsv', sep = '\t', header = T, row.names = 1)
metaZF <- read.table('zf.clusters.tsv', sep = '\t', header = T)

readnumbZF <- read.table('zf.snp_containing_reads_per_barcode', sep = '',col.names=c("SNPreads","barcode"))
separate(readnumbZF, barcode, sep = "-", c("barcodes", "originAgain")) -> readnumbZF

table(metaZF$status)
#doublet    singlet unassigned 
#334      10085        342 

table(metaZF$assignment)
#0  0/1  0/2    1  1/0  1/2    2  2/0  2/1 
#3286   85  126 1875   66   46 5068   94  115 

#from web summaries cell numbers for 
#D_1 = 2,193
#L_1 = 3,628
#M_1 = 6,016

#will keep unassinged 
#the cell origin is sorted in the barcode info barcode-# with the # being the origin
#if it has a D it's a doublet and if it has S its a doublet but from the same animal
separate(metaZF, barcode, sep = "-", c("barcodes", "origin")) -> metaZF
metaZF <- merge(metaZF,readnumbZF, by = "barcodes")

#going to clean things up a bit
metaZF -> cleaned.metaZF

#clean it up a bit by changing all #D in origin to doublet and all #/# in soup assignment to doublet
str_replace_all(cleaned.metaZF$assignment,"[:digit:]/[:digit:]", "Doublet") -> cleaned.metaZF$assignment
str_replace_all(cleaned.metaZF$origin,"[:digit:]D", "Doublet") -> cleaned.metaZF$origin

#remove lines with homotypic doublets since this is not a useful test for souporcell doublet detection
#troublet only detects doublets based on mixed genotypes
#cleaned.meta[cleaned.meta$origin != "[:digit:]S",] -> cleaned.metaTest
cleaned.metaZF[cleaned.metaZF$origin != "3S",] -> cleaned.metaZF
cleaned.metaZF[cleaned.metaZF$origin != "2S",] -> cleaned.metaZF
cleaned.metaZF[cleaned.metaZF$origin != "1S",] -> cleaned.metaZF

#here we're replacing the souporcell assignments with values matching the animal origin assigments
#that were found to correlate with them via an upset plot.
#str_replace_all(cleaned.meta$assignment,"[2]", "Animal#3") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[0]", "Animal#2") -> cleaned.meta$assignment
#str_replace_all(cleaned.meta$assignment,"[1]", "Animal#1") -> cleaned.meta$assignment
#str_replace_all(cleaned.metaZF$assignment,"[2]", "2") -> cleaned.metaZF$assignment
str_replace_all(cleaned.metaZF$assignment,"[0]", "3") -> cleaned.metaZF$assignment
#str_replace_all(cleaned.metaZF$assignment,"[1]", "1") -> cleaned.metaZF$assignment


#now combine souporcell and ground truth to generate bar graph for an upset plot
cleaned.metaZF$combined_assignment <-paste(cleaned.metaZF$origin,cleaned.metaZF$assignment, sep = '_')

#upset plot
ggplot(cleaned.metaZF, aes(x=reorder(combined_assignment, combined_assignment, function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme(axis.text.x = element_text(angle = 90))->upsetPlotZF

pdf("/Users/joecardiello/Documents/LeighLab/Data/220809_DemuxFiguresV2/Fig6_AxoZF/220827_ZF_Figure_Upset.pdf") 
#or do same as above but png, takes up way less room
print(upsetPlotZF)
dev.off()


#save for use in python
cleaned.metaZF[,c("status","assignment","combined_assignment", "origin","SNPreads")] -> df_ZF

names(df_ZF) <- c('soup_status', 'soup_assign', 'combined', 'origin',"SNPreads")

write.table(df_ZF, "/Users/joecardiello/Documents/LeighLab/Data/220827.ZF.ground.truth.soup.comparison.txt", sep = '\t', col.names = T)
```
