
library(GOstats)
library(DGCA) #wrapper for GOstats and visualization
library(HGNChelper)
library(org.Hs.eg.db)
library(ggplot2)
library(DBI)
options(stringsAsFactors = FALSE)

#make up a list of genes
data(darmanis) #scRNAseq data set of brain cell types, included in DGCA as example data
module_genes = list(
  mod1 = rownames(darmanis)[1:100],
  mod2 = rownames(darmanis)[90:190],
  mod3 = rownames(darmanis)[190:290],
  mod4 = rownames(darmanis)[330:360])
modules = stack(module_genes)
modules$ind = as.character(modules$ind)
str(modules)
head(modules)

#using your own gene list
gene_list = read.table("/Users/amckenz/Desktop/test_genes.txt")
module_genes = list(our_genes = gene_list$V1)
modules = stack(module_genes)
modules$ind = as.character(modules$ind)
str(modules)
head(modules)

#do GO analysis using GOstats
moduleGO_res = moduleGO(genes = modules$values, labels = modules$ind,
  universe = rownames(darmanis), pval_GO_cutoff = 1)
moduleGO_df = extractModuleGO(moduleGO_res)

#visualize the gene ontology categories
plotModuleGO(moduleGO_df, nTerms = 4, text_size = 8, coord_flip = TRUE)
