
library(rentrez)
options(stringsAsFactors = FALSE)
#library(bayesbio) #this is implemented as a function in library(bayesbio)

#example gene set
rowTerms = c("TP53", "TNF", "APOE", "EGFR", "ESR1", "IL6", "VEGFA",
  "MTHFR", "TGFB1", "IL10", "ACE", "BRCA1")
# rowTerms = c("CDC5L", "RELN", "TUBA1A", "CDC5L", "APOE", "ESR1")
colTerms = c("methylation", "immunity", "breast cancer", "Alzheimer")

disease_gene_mentions = data.frame(matrix(0, nrow = length(rowTerms),
   ncol = length(colTerms) + 1))

sleepTime = 0.01

for(i in 1:length(colTerms)){
	for(j in 1:length(rowTerms)){
    query = paste(colTerms[i], "AND", rowTerms[j], sep = " ")
    res = entrez_search(db="pubmed", term = query)
    disease_gene_mentions[j, i] = as.numeric(res$count)
    Sys.sleep(sleepTime)
	}
}

total_res = numeric(length(rowTerms))
for(j in 1:length(rowTerms)){
  res = entrez_search(db="pubmed", term = rowTerms[j])
  total_res[j] = as.numeric(res$count)
  Sys.sleep(sleepTime)
}

rownames(disease_gene_mentions) = rowTerms
disease_gene_mentions[, length(colTerms) + 1] = total_res
colnames(disease_gene_mentions) = c(colTerms, "Total_Mentions")
disease_gene_mentions = disease_gene_mentions[rev(order(disease_gene_mentions[ , length(colTerms) + 1])), ]
# disease_gene_mentions = disease_gene_mentions[rev(order(disease_gene_mentions[ , "Alzheimer"])), ]

########
#creating a year-by-year plot of a particular term

#code from https://ropensci.org/tutorials/rentrez_tutorial.html
search_year <- function(year, term){
    query = paste(term, "AND (", year, "[PDAT])")
    entrez_search(db="pubmed", term=query, retmax=0)$count
}

year = 1970:2016
search_term = "BRCA1"
term_results = sapply(year, search_year, term=search_term, USE.NAMES=FALSE)
paper_total = sapply(year, search_year, term="", USE.NAMES=FALSE)
term_relative = term_results/paper_total
term_relative = term_relative/sum(term_relative)

plot(year, term_relative, type='b', ylab = "Proportion of Searches",
  main=paste0("Proportion of PubMed Searches for ", search_term))
