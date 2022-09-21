# SNPdemuxPaper
for depositing scripts used in benchmarking SNP demux data

With 10X Cellranger count and multi generated files, as well as souporcell generated clusters.tsv files, we first read these into R for filterring,
correlation analysis and subsequent assignment label renaming, UMAP plotting, and assembling/saving of dataframes for python analysis.

R scripts for this analysis are in this github page and are labelled by figure they were used for, and species. These R scripts were used to produce the 
dataframe files needed for python analyses posted bellow. 

Python analysis took place in individual google colab notebooks here:
Zebrafish, axolotl, and green monkey: https://colab.research.google.com/drive/1yXzE3WJ05hEJKdy7owiXOpCjYvUm4DDJ#scrollTo=xoYMouNlv1JK
Xenopus: https://colab.research.google.com/drive/1lO4ny8Uv9n1lPIbHmZKFPyxWI7gjmZPs?usp=sharing#scrollTo=5NbAyigRiBsx
Pleurodeles: https://colab.research.google.com/drive/1Zbbpi1WwfKGwrrFhuecrHE3lSjP0Pzsz#scrollTo=dJTlyRO8ipjf
Pleurodeles and Notophthalmus cellplex: https://colab.research.google.com/drive/12ZNvvfiUt3DL6UpTg8BjQMSd4Y-6yM3u#scrollTo=o-G_FjoCYO_f
Pleurodeles and Notophthalmus barnyard:https://colab.research.google.com/drive/1JS8kRUGAioDM2IvYa4oBDzNMKLxBhHV_#scrollTo=dJTlyRO8ipjf
