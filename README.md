# SNPdemuxPaper
Scripts used in Cardiello et al. 2022 Biorxiv "Accurate genotype-based demultiplexing of single cell RNA sequencing samples from non-human animals" 

With 10X Cellranger count and multi generated files, as well as souporcell generated clusters.tsv files, we first read these into R for filtering,
correlation analysis and subsequent assignment label renaming, UMAP plotting, and assembling/saving of dataframes for python analysis.

R scripts for this analysis are in this github page and are labelled by figure and species they were used. These R scripts were used to produce the 
dataframe files needed for python analyses posted below. 

Python analysis took place in individual google colab notebooks here:

Zebrafish, axolotl, and green monkey (Figure 1, and Supplemental Figure 1): https://colab.research.google.com/drive/1yXzE3WJ05hEJKdy7owiXOpCjYvUm4DDJ#scrollTo=xoYMouNlv1JK

Xenopus (Figure 2): https://colab.research.google.com/drive/1lO4ny8Uv9n1lPIbHmZKFPyxWI7gjmZPs?usp=sharing#scrollTo=5NbAyigRiBsx

Pleurodeles (Figure 3): https://colab.research.google.com/drive/1Zbbpi1WwfKGwrrFhuecrHE3lSjP0Pzsz#scrollTo=dJTlyRO8ipjf
We have also included the file plasmids.add.to.ref.fa which can be appended to a Pleurodeles reference to map the fluorscent genes for each of the fluorescent Pleurodeles strains.

Pleurodeles and Notophthalmus cellplex (Figure 4): https://colab.research.google.com/drive/12ZNvvfiUt3DL6UpTg8BjQMSd4Y-6yM3u#scrollTo=o-G_FjoCYO_f

Pleurodeles and Notophthalmus barnyard (Supplemental Figure 4): https://colab.research.google.com/drive/1JS8kRUGAioDM2IvYa4oBDzNMKLxBhHV_#scrollTo=dJTlyRO8ipjf
