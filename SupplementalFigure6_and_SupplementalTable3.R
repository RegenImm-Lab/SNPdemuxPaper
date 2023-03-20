library(tidyverse)
library(ggplot2)
library(tidyr)
library(Seurat)
library(patchwork)

#analysis of pool of 30 zebrafish embryos, overlaying souporcell assignments on clustering results

#load in souporcell outputs
soup.data <- read.csv(file = 'clusters.tsv', header = T, sep = '\t', row.names = 1)

#subset to make a table of the singlet assignments
soup.data[soup.data$status == 'singlet',] -> singlets
table(singlets$assignment)

#now into Seurat, need to specify this or have setwd to location of data
zf30.data <- Read10X(data.dir = 'filtered_feature_bc_matrix/')

#per paper (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0254024#sec002): The gene counts matrix was loaded into Seurat and a Seurat object was created by filtering cells 
#which only expressed more than 200 genes and filtering genes that were expressed in at least 3 cells
zf30 <- CreateSeuratObject(counts = zf30.data, project = "zf30", min.cells = 3, min.features = 200, meta.data = soup.data)
zf30
#20986 features across 25240 samples within 1 assay 

#will follow paper analysis
#Additionally, as an extra quality-control step, cells were filtered out (excluded) based on the following 
#criteria: <400 or >2,500 unique genes expressed, or >5% of counts mapping to the mitochondrial genome. 
#This resulted in 20,279 cells in the dataset.

#find percent mito
zf30[["percent.mt"]] <- PercentageFeatureSet(zf30, pattern = "^mt-")

#quick inspection
VlnPlot(zf30, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(zf30, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zf30, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#paper filters quite harshly and was old version of cellranger so likely 
#intron counts are factoring in to some higher mitochondrial and nFeatures here
#will instead do 10% mito and 3000 features
zf30 <- subset(zf30, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)
zf30
#20986 features across 19665 samples within 1 assay 
#Active assay: RNA (20986 features, 0 variable features)
#this diverges from paper, a bit but if we keep the same cutoffs then its about 3k cells less
#whereas this is <1k different

#normalize
zf30 <- NormalizeData(zf30, normalization.method = "LogNormalize", scale.factor = 10000)

zf30 <- FindVariableFeatures(zf30, selection.method = "vst", nfeatures = 2000)

#10 most hvg, standard seurat pipeline
top10 <- head(VariableFeatures(zf30), 10)

# plot variable genes with + without labels
plot1 <- VariableFeaturePlot(zf30)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(zf30)
zf30 <- ScaleData(zf30, features = all.genes)

zf30 <- RunPCA(zf30, features = VariableFeatures(object = zf30))

#paper used elbow so we will see that too
ElbowPlot(zf30, ndims = 50)

#they took 40 PCs, found 22 clusters using resolution of 0.5 (louvian based)

zf30 <- FindNeighbors(zf30, dims = 1:40)
zf30 <- FindClusters(zf30, resolution = 0.5)
#so we get 21 clusters and they had 22, pretty close considering differences in alignment and filtering

#run umap
zf30 <- RunUMAP(zf30, dims = 1:40)

#plot umap
DimPlot(zf30, reduction = "umap")
zf30$assignment
#look at doublets...still a lot in here which can be removed given we have this info
DimPlot(zf30, group.by = 'status')

#are there clusters derived from one of few individuals? 

subset(x = zf30, subset = status == "singlet") -> zf.30.singlets

DimPlot(zf.30.singlets, group.by = 'assignment')

data.frame(rbind(table(zf.30.singlets$seurat_clusters, zf.30.singlets$assignment))) -> cluster.composition
#supplemental table 3, this was then opened in Excel and conditionally formatted 
write.table(cluster.composition, 'cluster.composition.txt', quote = F, sep = '\t', row.names = T, col.names = T)
cluster.composition

#how many clusters have at least 5 cells from each individual?
#how many rows (so clusters) have less than 30 biological replicates
cluster.composition$under.five.rows <- apply(cluster.composition, 1, function(x) length(which(x >= 5)))
#[1] 30 30 30 30 30 30 25 29 27 27 27 24 18 19 14 14 12 11 11  8  7
#conclusion: While 6 clusters have biological replicates from 30 individuals, the remaining 15 do not have more than 5 cells
#lowest is 7, meaning that this cluster only has more than 5 cells from 7 replicates

#what about a more stringent zero?
cluster.composition$zero.rows <- apply(cluster.composition, 1, function(x) length(which(x=="0")))
#output  [1] 0 0 0 0 0 0 0 0 1 1 1 0 1 2 5 4 6 4 6 9 4
#9/21 have a zero which means 9/21 have at least 1 cell from each zf, the lowest is 9, meaning that
#21 of the 30 biological replicates are in this sample

#let's look at doublets/cluster
data.frame(rbind(table(zf30$seurat_clusters, zf30$status))) -> doublet.cluster.composition
doublet.cluster.composition

#let's look at % doublets/population
doublet.cluster.composition$singlet + doublet.cluster.composition$doublet -> doublet.cluster.composition$total
(doublet.cluster.composition$doublet/doublet.cluster.composition$total)*100 -> doublet.cluster.composition$percent.doublet
doublet.cluster.composition$percent.doublet
#[1] 12.9978388 14.4786976 25.0217960 23.3752621 15.0031786 18.5981308 29.8524404 15.1915456 27.1056662 17.0103093
#[11]  0.7692308 17.1366594 21.0937500 28.9772727 28.7128713 20.5980066 19.0871369 30.2752294 22.7272727 27.4509804
#[21] 10.6796117
#conclusion: all clusters have some doublets intermixed, this ranges from very low 0.8% doublets to 30%, which means that 
#one out of every 3 cells in that cluster is a doublet

doublets <- as.numeric(doublet.cluster.composition$percent.doublet)
summary(doublets)

#supplemental Figure 6B
pdf(file = 'doublet.composition.pdf', width = 4, height = 4)
#graph this
boxplot(doublet.cluster.composition$percent.doublet)
stripchart(doublet.cluster.composition$percent.doublet,             
           method = "jitter", 
           pch = 19,          
           col = 4,           
           vertical = TRUE,  
           add = TRUE) 
dev.off()

clusters <- as.numeric(table(singlets$assignment))
summary(clusters)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#275.0   405.5   532.5   644.6   853.5  1517.0 

#supplemental Figure 6A
pdf(file = 'cells.per.animal.pdf', width = 4, height = 4)
#quick qc graph of the contribution/animal
boxplot(clusters)
stripchart(clusters,             
           method = "jitter", 
           pch = 19,          
           col = 4,           
           vertical = TRUE,  
           add = TRUE)                            

dev.off()

sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin20.6.0 (64-bit)
Running under: macOS Big Sur 11.7.4

Matrix products: default
LAPACK: /usr/local/Cellar/r/4.2.2/lib/R/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] patchwork_1.1.2    forcats_0.5.2      purrr_1.0.1        readr_2.1.3        tibble_3.1.8      
[6] tidyverse_1.3.2    tidyr_1.2.1        Matrix_1.5-3       ggplot2_3.4.0      stringr_1.5.0     
[11] dplyr_1.0.10       SeuratObject_4.1.3 Seurat_4.3.0      

loaded via a namespace (and not attached):
  [1] googledrive_2.0.0      Rtsne_0.16             ggbeeswarm_0.7.1       colorspace_2.0-3      
[5] deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4         fs_1.5.2              
[9] rstudioapi_0.14        spatstat.data_3.0-0    leiden_0.4.3           listenv_0.9.0         
[13] farver_2.1.1           ggrepel_0.9.2          lubridate_1.9.0        fansi_1.0.3           
[17] xml2_1.3.3             codetools_0.2-18       splines_4.2.2          polyclip_1.10-4       
[21] jsonlite_1.8.4         broom_1.0.2            ica_1.0-3              dbplyr_2.2.1          
[25] cluster_2.1.4          png_0.1-8              uwot_0.1.14            shiny_1.7.4           
[29] sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.2         httr_1.4.4            
[33] backports_1.4.1        assertthat_0.2.1       fastmap_1.1.0          lazyeval_0.2.2        
[37] gargle_1.2.1           cli_3.6.0              later_1.3.0            htmltools_0.5.4       
[41] tools_4.2.2            igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
[45] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.9             scattermore_0.8       
[49] cellranger_1.1.0       vctrs_0.5.1            spatstat.explore_3.0-5 nlme_3.1-161          
[53] progressr_0.13.0       lmtest_0.9-40          spatstat.random_3.0-1  globals_0.16.2        
[57] rvest_1.0.3            timechange_0.2.0       mime_0.12              miniUI_0.1.1.1        
[61] lifecycle_1.0.3        irlba_2.3.5.1          googlesheets4_1.0.1    goftest_1.2-3         
[65] future_1.30.0          MASS_7.3-58.1          zoo_1.8-11             scales_1.2.1          
[69] hms_1.1.2              promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.2        
[73] RColorBrewer_1.1-3     reticulate_1.27        pbapply_1.6-0          gridExtra_2.3         
[77] ggrastr_1.0.1          stringi_1.7.12         rlang_1.0.6            pkgconfig_2.0.3       
[81] matrixStats_0.63.0     lattice_0.20-45        ROCR_1.0-11            tensor_1.5            
[85] htmlwidgets_1.6.1      labeling_0.4.2         cowplot_1.1.1          tidyselect_1.2.0      
[89] parallelly_1.33.0      RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
[93] R6_2.5.1               generics_0.1.3         DBI_1.1.3              haven_2.5.1           
[97] pillar_1.8.1           withr_2.5.0            fitdistrplus_1.1-8     survival_3.5-0        
[101] abind_1.4-5            sp_1.5-1               future.apply_1.10.0    modelr_0.1.10         
[105] crayon_1.5.2           KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3   
[109] plotly_4.10.1          tzdb_0.3.0             readxl_1.4.1           grid_4.2.2            
[113] data.table_1.14.6      reprex_2.0.2           digest_0.6.31          xtable_1.8-4          
[117] httpuv_1.6.8           munsell_0.5.0          beeswarm_0.4.0         viridisLite_0.4.1     
[121] vipor_0.4.5   