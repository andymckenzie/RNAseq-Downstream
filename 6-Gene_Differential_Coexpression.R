#some taken from the DGCA vignette http://htmlpreview.github.io/?https://github.com/andymckenzie/DGCA/blob/master/inst/doc/DGCA.html

library(DGCA, quietly = TRUE)

options(stringsAsFactors = FALSE)

data(darmanis)
data(design_mat)

str(darmanis)
darmanis_mean_filtered = filterGenes(darmanis,
  filterTypes = "central", filterCentralType = "mean", filterCentralPercentile = 0.75)
str(darmanis_mean_filtered)

#Calculate the gene-level differential correlations on the filtered matrix
ddcor_res = ddcorAll(inputMat = darmanis_mean_filtered, design = design_mat,
  compare = c("oligodendrocyte", "neuron"), corrType = "spearman",
  adjust = "perm", heatmapPlot = FALSE, nPerm = 10)
str(ddcor_res)
head(ddcor_res)

#Visualization of gene pairs that are differential correlated between conditions
library(ggplot2, quietly = TRUE)
#remove one outlier sample before visualization
darmanis = darmanis[ , -which.max(darmanis["COX6A1", ])]
design_mat = design_mat[-which.max(darmanis["COX6A1", ]), ]
plotCors(inputMat = darmanis, design = design_mat, compare = c("oligodendrocyte", "neuron"), geneA = "RTN4", geneB = "COX6A1")
plotCors(inputMat = darmanis, design = design_mat, compare = c("oligodendrocyte", "neuron"), geneA = "CALM2", geneB = "UBB")

#Visualizing the overall heatmap of correlations in each condition
library(gplots, quietly = TRUE)
darmanis_top =  filterGenes(darmanis,
  filterTypes = c("central", "dispersion"), filterCentralPercentile = 0.75,
  filterDispersionPercentile = 0.75)
ddcor_res = ddcorAll(inputMat = darmanis_top, design = design_mat,
  compare = c("oligodendrocyte", "neuron"),
  adjust = "none", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")

#Gene-trait differential correlation analysis
data(ages_darmanis)
rownames_darmanis = rownames(darmanis)
darmanis_with_traits = rbind(ages_darmanis, darmanis)
rownames(darmanis_with_traits) = c("age", rownames_darmanis)
ddcor_res = ddcorAll(inputMat = darmanis_with_traits, design = design_mat,
  compare = c("oligodendrocyte", "neuron"), corrType = "spearman",
  adjust = "perm", nPerm = 10, splitSet = "age")
head(ddcor_res)
#plot some of the differential correlations here 
