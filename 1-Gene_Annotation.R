
library(rentrez)
options(stringsAsFactors = FALSE)

#example gene set
rowTerms = c("APOE", "VHL", "EGFR")
# rowTerms = c("CDC5L", "RELN", "TUBA1A", "CDC5L", "APOE", "ESR1")
colTerms = c("cancer", "alzheimer", "multiple sclerosis")

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

year = c(1980:2016)
term = "MTOR"

term_results = vector()

for(j in 1:length(year)){
#   query = paste(term, "AND (", year, "[PDAT])")
#   res = entrez_search(db="pubmed", term = rowTerms[j])
  term_results[j] = as.numeric(search_year(year[j], term))
  Sys.sleep(sleepTime)
}

# term_results = sapply(year, search_year, term=search_term, USE.NAMES=FALSE)
paper_total = sapply(year, search_year, term="", USE.NAMES=FALSE)
term_relative = term_results/paper_total
# term_relative = term_relative/sum(term_relative)

plot(year, term_relative, type='b', ylab = "Proportion of Searches",
  main=paste0("Proportion of PubMed Searches for ", term))

plot(year, term_relative, type='b', ylab = "Proportion of Searches",
  main=paste0("Proportion of PubMed Searches for ", term), xaxt="n")
axis(1, at = 1980:2016, las=2)
