#install required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSEABase")
BiocManager::install('affy')

#calling required libraries

BiocManager::install("hgu133plus2.db")
library(BiocManager)
library(hgu133plus2.db)
library(tidyverse)
library(affy)
library(GSEABase)

#set working directories
setwd('/project/bf528/project_1/data/')

#read your csv file
biological_results <- read.csv ('differential_expression_results.csv', col.names = c('PROBEID', 't', 'p', 'p_adjust'))

#Using the select() function of the bioconductor package hgu133plus2.db, map the probeset IDs to gene symbols by specifying the appropriate key and column arguments. Some probeset IDs map to the same gene symbol, so reason about and pick a rationale for choosing which probeset ID to use as representative. Add an additional column to the differential expression results that contains one symbol for each probeset ID.
genesymbol <- AnnotationDbi::select(hgu133plus2.db, biological_results$PROBEID, c('SYMBOL'))
head(biological_results)
head(genesymbol)
View(genesymbol)

#collapser is a function that combines genesymbols that map to 2 different probeids using '|' 
collapser <- function(a) {a %>% unique %>% sort %>% paste(collapse = '|') }
#using collapser, we combine multiple genesymbols with same probeid
genesymbol <- genesymbol %>% group_by(PROBEID) %>% summarise_each(funs(collapser)) %>% ungroup
View(genesymbol)

#storing and merging biological results with genesymbol 
merge_symbol <- merge(biological_results, genesymbol, on = 'PROBEID')
View(merge_symbol)

#removes blank rows
merge_symbol <- merge_symbol[!(is.na(merge_symbol$SYMBOL) | merge_symbol$SYMBOL == ""), ]
View(merge_symbol)

#to find no. of ids with same genesymbol
countsignificant <- merge_symbol %>% group_by(SYMBOL) %>% count() %>% filter(n>=2)
countsignificant[order(countsignificant$n), ]
View(countsignificant)

#finding the most significant with common genesymbols
idiff <- data.frame(PROBEID = character(), t = numeric(), p = numeric(), p_adjust = numeric(), SYMBOL = character() )
for (i in countsignificant$SYMBOL) {
  x <- merge_symbol[merge_symbol$SYMBOL==i, ]
  x <- x[x$p_adjust == min(x$p_adjust), ]
  merge_symbol <- merge_symbol[!merge_symbol$SYMBOL == i,]
  idiff <- rbind(idiff, x)
}
merge_symbol <- rbind(merge_symbol, idiff)

#load gene sets 
hallmarks <- getGmt('h.all.v7.2.symbols.gmt')
GO <- getGmt('c5.go.v7.2.symbols.gmt')
KEGG <- getGmt('c2.cp.kegg.v7.2.symbols.gmt')

#get geneset length
length(hallmarks)      #50
length(GO)             #10271    
length(KEGG)           #186

#order the values in decreasing values of t-stats 
merge_symbol <- merge_symbol[order(merge_symbol$t, decreasing = TRUE), ]
head(merge_symbol)
tail(merge_symbol)
View(merge_symbol)

#top 1000 up regulated and down regulated
up_1000 <- head(merge_symbol, n = 1000)
down_1000 <- tail(merge_symbol, n = 1000)
View(up_1000)

#selecting top 10 for report
up_10 <- head(up_1000, n = 10)
down_10 <- tail(down_1000, n = 10)

#store the top10 up and down regulated genes
write.csv(up_10, "10_upregulated_genes.csv")
write.csv(down_10, "10_downregulated_genes.csv")

#extract the genes that were not expressed
not_diffexp_up <- subset(merge_symbol, !merge_symbol$SYMBOL %in% up_1000$SYMBOL)
not_diffexp_down <- subset(merge_symbol, !merge_symbol$SYMBOL %in% down_1000$SYMBOL)

#define function to create contigency table
fishertest <- function(gl, gs, nde)           #gl = genelist, gs= geneset, nde= not differentially expressed
{ diffexp_ings <- length(intersect(gl,gs))    #diffexp_ings = differentially expressed genes that are present in geneset
diffexp_notgs <- length(gl) - diffexp_ings    #diffexp_notgs = differentially expressed genes that are not present in geneset 
notde_ings <- length(intersect(nde,gs))       #notde_ings = not expressed but present in geneset
notde_notgs <- length(nde) - notde_ings       #notde_notgs = not differentially expressed and not in geneset
return(c(diffexp_ings,diffexp_notgs,notde_ings,notde_notgs))}   #returns the fishertest values

#stores results of fisher test for hallmark geneset 
hallmarks_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

#stores the results for hallmark geneset comparison in separate data frame using for loop
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}
View(hallmarks_results)
hallmarks_results <- hallmarks_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(hallmarks_results)

#stores results of fisher test for kegg geneset
kegg_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

##stores the results for kegg geneset comparison in separate data frame using for loop
for (i in 1:length(KEGG))
{
  geneid <- geneIds(KEGG[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}

kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(kegg_results)    #to view result table

#stores results of fisher test for kegg geneset
go_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

##stores the results for go geneset comparison in separate data frame using for loop
for (i in 1:length(GO))
{
  geneid <- geneIds(GO[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  go_results[nrow(go_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  go_results[nrow(go_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}

go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(go_results)   #to view result table

#adjusting the pvalue usinf benjamini hochberg method and storing the data in seperate file
go_results$BH <- p.adjust(go_results$pvalue, method = "BH", n = length(go_results$pvalue))
write.csv(go_results, "final_go.csv")     

kegg_results$BH <- p.adjust(kegg_results$pvalue, method = "BH", n = length(kegg_results$pvalue))
write.csv(kegg_results, "final_kegg.csv")

hallmarks_results$BH <- p.adjust(hallmarks_results$pvalue, method = "BH", n = length(hallmarks_results$pvalue))
write.csv(hallmarks_results, "final_hallmarks.csv")