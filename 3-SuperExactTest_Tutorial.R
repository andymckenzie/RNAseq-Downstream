#install.packages("SuperExactTest")
library(SuperExactTest)
options(stringsAsFactors = FALSE)

###############
# Test data using letters to show how intersections work

testA = letters[1:10]
testB = letters[7:20]
testC = letters[10:26]

meta_list = list(testA, testB, testC)

names(meta_list) = c("testA", "testB", "testC")

meta_result = supertest(meta_list, n = 26)
summary(meta_result)

#examine the result of using different "universe sizes"
meta_result = supertest(meta_list, n = 100)
summary(meta_result)

meta_result = supertest(meta_list, n = 1000)
summary(meta_result)

###########
# Compare differential expression from two public data sets using GEO2R

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE84495", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL11533", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00001111XXXX22223333XXXX444445555XXXXXXXXXXXXXXXXX",
        "XX1X223XX45XXXX1X23X445XXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)

#NB: added some of the below in manually
colnames(design) <- c("Control_6MO", "CR_6MO", "Control_12MO", "CR_12MO", "Control_24MO", "CR_24MO")
fit <- lmFit(gset, design)

significance_threshold = 0.01

cont.matrix_6mo <- makeContrasts(Control_6MO-CR_6MO, levels=design)
fit2_6mo <- contrasts.fit(fit, cont.matrix_6mo)
fit2_6mo <- eBayes(fit2_6mo, 0.01)
tT_6mo <- topTable(fit2_6mo, adjust="fdr", sort.by="B", number=nrow(ex))
tT_6mo <- tT_6mo[ , c("ID", "Gene.title", "Gene.symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
tT_6mo_sig <- tT_6mo[tT_6mo$adj.P.Val < significance_threshold, ]

cont.matrix_12mo <- makeContrasts(Control_12MO-CR_12MO, levels=design)
fit2_12mo <- contrasts.fit(fit, cont.matrix_12mo)
fit2_12mo <- eBayes(fit2_12mo, 0.01)
tT_12mo <- topTable(fit2_12mo, adjust="fdr", sort.by="B", number=nrow(ex))
tT_12mo <- tT_12mo[ , c("ID", "Gene.title", "Gene.symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
tT_12mo_sig <- tT_12mo[tT_12mo$adj.P.Val < significance_threshold, ]

cont.matrix_24mo <- makeContrasts(Control_24MO-CR_24MO, levels=design)
fit2_24mo <- contrasts.fit(fit, cont.matrix_24mo)
fit2_24mo <- eBayes(fit2_24mo, 0.01)
tT_24mo <- topTable(fit2_24mo, adjust="fdr", sort.by="B", number=nrow(ex))
tT_24mo <- tT_24mo[ , c("ID", "Gene.title", "Gene.symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
tT_24mo_sig <- tT_24mo[tT_24mo$adj.P.Val < significance_threshold, ]

cr_meta_list = list(unique(tT_6mo_sig$Gene.symbol), unique(tT_12mo_sig$Gene.symbol),
  unique(tT_24mo_sig$Gene.symbol))
names(cr_meta_list) = c("6MO_Control_vs_CR", "12MO_Control_vs_CR", "24MO_Control_vs_CR")
cr_meta_result = supertest(cr_meta_list, n = length(unique(tT_6mo$Gene.symbol)))
summary(cr_meta_result)

##########
# If you have time, try using GEO2R on a data set(s)/question that you are interested in!
