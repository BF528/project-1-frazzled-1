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
#setwd('/project/bf528/project_1/data/')
getwd()
#read your csv file
biological_results <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv", row.names = 1, header = TRUE)
biological_results <- biological_results %>% arrange(desc(t))
view(biological_results)

#Using the select() function of the bioconductor package hgu133plus2.db, map the probeset IDs to gene symbols by specifying the appropriate key and column arguments.
#Add an additional column to the differential expression results that contains one symbol for each probeset ID.
de_matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(row.names(biological_results)), columns = ("SYMBOL"))
# Remove duplicated matches
de_duplicatedmatches <- de_matches[!duplicated(de_matches[1]),]
# Combine symbols with initial de_data to add symbol column
biological_results <- cbind(de_duplicatedmatches, biological_results)
# Keeping only the probes with significant padj (adjust p-values)
biological_results <- biological_results %>%
  group_by(SYMBOL) %>%
  filter(padj == min(padj)) %>%
  ungroup(SYMBOL)

#load gene sets 
hallmarks <- getGmt('h.all.v7.5.1.symbols.gmt')
GO <- getGmt('c5.go.v7.5.1.symbols.gmt')
KEGG <- getGmt('c2.cp.kegg.v7.5.1.symbols.gmt')

#get geneset length
length(hallmarks)      #50
length(GO)             #10271    
length(KEGG)           #186

# Number of gene sets per collection
cat("Number of gene sets in Hallmark: ", length(names(hallmarks)))
cat("Number of gene sets in GO: ", length(names(GO)))
cat("Number of gene sets in KEGG: ", length(names(KEGG)))

# Remove missing values
biological_results <- biological_results[!is.na(biological_results$SYMBOL),]

# Gather top 1000 up-regulated and down-regulated genes 
upreg_1000 <- head(biological_results, 1000)
downreg_1000 <- tail(biological_results, 1000)

# Gather top 10 up-regulated and down-regulated genes from top 1000 list
upreg_10 <- head(upreg_1000, 10)
downreg_10 <- tail(downreg_1000, 10)

#store the top10 up and down regulated genes
write.csv(upreg_10, "10_upregulated_genes.csv")
write.csv(downreg_10, "10_downregulated_genes.csv")

#extract the genes that were not expressed
not_diffexp_up <- subset(biological_results, ! biological_results $SYMBOL %in% upreg_1000 $ SYMBOL)
not_diffexp_down <- subset(biological_results, ! biological_results $SYMBOL %in% downreg_1000 $ SYMBOL)

#define function to create contigency table
fishertest <- function(gl, gs, nde)           #gl = genelist, gs= geneset, nde= not differentially expressed
{ 
  diffexp_ings <- length(intersect(gl,gs))    #diffexp_ings = differentially expressed genes that are present in geneset
diffexp_notgs <- length(gl) - diffexp_ings    #diffexp_notgs = differentially expressed genes that are not present in geneset 
notde_ings <- length(intersect(nde,gs))       #notde_ings = not expressed but present in geneset
notde_notgs <- length(nde) - notde_ings       #notde_notgs = not differentially expressed and not in geneset
return(c(diffexp_ings,diffexp_notgs,notde_ings,notde_notgs))  #returns the fishertest values
}  

#stores results of fisher test for hallmark geneset 
hallmarks_results <- data.frame(setname = character(), 
                                pvalue = numeric(), estimate = numeric(), 
                                exp = character(), stringsAsFactors = FALSE)

#stores the results for hallmark geneset comparison in separate data frame using for loop
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(upreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')}
#View(hallmarks_results)
hallmarks_results <- hallmarks_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))
View(hallmarks_results)

#stores results of fisher test for kegg geneset
kegg_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

##stores the results for kegg geneset comparison in separate data frame using for loop
for (i in 1:length(KEGG))
{
  geneid <- geneIds(KEGG[i])
  fisher_up <- fishertest(upreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
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
  fisher_up <- fishertest(upreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
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

# Statistically significant and enriched genesets
# Kegg
kegg_sig <- kegg_results[kegg_results$pvalue<0.05,]
kegg_enr_num <- length(kegg_sig$gene_set)
cat("Statistically enriched genesets in KEGG: ", kegg_enr_num, "\n")

# GO:
go_sig <- go_results[go_results$pvalue<0.05,]
go_enr_num <- length(go_sig$gene_set)
cat("Statistically enriched genesets in GO: ", go_enr_num, "\n")

# Hallmark:
hm_sig <- hallmarks_results[hallmarks_results$pvalue<0.05,]
hm_enr_num <- length(hm_sig$gene_set)
cat("Statistically enriched genesets in Hallmark: ", hm_enr_num, "\n")

# Top 3 genesets for each:
top3_kegg <- slice_min(kegg_results, order_by = pvalue, n=3)
top3_go <- slice_min(go_results, order_by = pvalue, n=3)
top3_hallmark <- slice_min(hallmarks_results, order_by = pvalue, n=3)
top_3 <- rbind(top3_kegg, top3_go, top3_hallmark)
print(top_3)
write.csv(top_3, file="geneset_top_results.csv")
