
library(WGCNA)
#WGCNA code adapted from https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
library(DGCA) #for example data
options(stringsAsFactors = FALSE)

#######

data(darmanis) #scRNAseq data set of brain cell types, included in DGCA as example data
gnxp = t(darmanis) #need to convert from columns-samples, genes-rows to vice versa

gsg = goodSamplesGenes(gnxp, verbose = 3);

#choose soft threshold and plot
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(gnxp, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
  xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 4

######
#one-step

net = blockwiseModules(gnxp, power = softPower,
                       TOMType = "unsigned", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
#
moduleLabels = net$colors;
moduleColors = labels2colors(moduleLabels)
modMembers = data.frame(Gene = colnames(gnxp), Module = moduleColors)

#######
# multi-step

adjacency = adjacency(gnxp, power = softPower)
str(adjacency)
TOM = TOMsimilarity(adjacency) #Topological Overlap Matrix
str(TOM)
dissTOM = 1-TOM
str(dissTOM)

geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)

minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors")

modMembersStep = data.frame(Gene = colnames(gnxp), Module = dynamicColors)

#######
# GO analysis

library(GOstats)
library(DGCA) #wrapper for GOstats and visualization
library(HGNChelper)
library(org.Hs.eg.db)
library(ggplot2)
library(DBI)

#do GO analysis using GOstats
moduleGO_res = moduleGO(genes = modMembersStep$Gene, labels = modMembersStep$Module,
  universe = colnames(gnxp), pval_GO_cutoff = 1)
moduleGO_df = extractModuleGO(moduleGO_res)
str(moduleGO_df)

#visualize the gene ontology categories
plotModuleGO(moduleGO_df, nTerms = 2, text_size = 5, coord_flip = FALSE)

modMembersStep[modMembersStep$Module == "yellow", ]

#############
# export to cytoscape

exportNetworkToCytoscape(adjacency, edgeFile = "edge_result.txt",
  nodeFile = "node_result.txt", weighted= TRUE, threshold = 0.1,
  nodeNames = rownames(adjacency))

write.table(modMembersStep, sep = "\t", file = "module_members.txt", row.names = F, quote = F) 

#########
# Part 2: Differential Coexpression of the networks between cell types
# need to leverage the design matrix for cell types

data(design_mat)
str(design_mat)

moduleDC_res = moduleDC(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"), genes = modMembersStep$Gene,
  labels = modMembersStep$Module, nPerm = 50, number_DC_genes = 3,
  dCorAvgMethod = "median")
head(moduleDC_res)
