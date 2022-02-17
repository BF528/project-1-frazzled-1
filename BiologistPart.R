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
biological_results <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv",
                               row.names = 1, header = TRUE)
biological_results <- biological_results %>% arrange(desc(t))
view(biological_results)

#Using the select() function of the bioconductor package hgu133plus2.db, 
#map the probeset IDs to gene symbols by specifying the appropriate key and column arguments.
#Add an additional column to the differential expression results that contains one symbol for each probeset ID.
de_matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character
                                    (row.names(biological_results)), 
                                    columns = ("SYMBOL"))
# Remove duplicated matches
de_duplicatedmatches <- de_matches[!duplicated(de_matches[1]),]

# Combine symbols with initial de_data to add symbol column
biological_results <- cbind(de_duplicatedmatches, biological_results)

# Keeping only the probes with significant padj (adjust p-values)
biological_results <- biological_results %>%
  group_by(SYMBOL) %>%
  filter(padj == min(padj)) %>%
  ungroup(SYMBOL)

#loading all the required gene sets 
hallmarks <- getGmt('h.all.v7.5.1.symbols.gmt')
GO <- getGmt('c5.go.v7.5.1.symbols.gmt')
KEGG <- getGmt('c2.cp.kegg.v7.5.1.symbols.gmt')

# Number of gene sets per collection
cat("Number of gene sets in Hallmark: ", length(names(hallmarks)))
cat("Number of gene sets in GO: ", length(names(GO)))
cat("Number of gene sets in KEGG: ", length(names(KEGG)))

#Remove all the missing values
biological_results <- biological_results[!is.na(biological_results $ SYMBOL),]

#Collect the top 1000 up-regulated and down-regulated genes 
upreg_1000 <- head(biological_results, 1000)
downreg_1000 <- tail(biological_results, 1000)

#Collect the top 10 up-regulated and down-regulated genes from above 1000 list
upreg_10 <- head(upreg_1000, 10)
downreg_10 <- tail(downreg_1000, 10)

#store the top10 up and down regulated genes
write.csv(upreg_10, "10_upregulated_genes.csv")
write.csv(downreg_10, "10_downregulated_genes.csv")

#extract the genes that were not expressed
not_diffexp_up <- subset(biological_results, ! biological_results $ SYMBOL 
                         %in% upreg_1000 $ SYMBOL)
not_diffexp_down <- subset(biological_results, ! biological_results $ SYMBOL
                           %in% downreg_1000 $ SYMBOL)

#defining a function to create contingency table
#gl = gene-list, gs= gene-set, nde= not differential expressed
#diffexp_in_gs = differentially expressed genes that are present in gene-set
#diffexp_not_gs = differentially expressed genes that are not present in gene-set
#notexp_in_gs = not expressed but present in gene-set
#notexp_not_gs = not differentially expressed and not in gene-set

fishertest <- function(gl, gs, nde)           
{ 
  diffexp_in_gs <- length(intersect(gl,gs))    
  diffexp_not_gs <- length(gl) - diffexp_in_gs   
  notexp_in_gs <- length(intersect(nde,gs))      
  notexp_not_gs <- length(nde) - notexp_in_gs      
  return(c(diffexp_in_gs,diffexp_not_gs,notexp_in_gs,notexp_not_gs))
}  

#stores results of fisher test for hallmark geneset 
hallmarks_results <- data.frame(setname = character(), 
                                pvalue = numeric(), estimate = numeric(), 
                                exp = character(), stringsAsFactors = FALSE)

#To store the results for hallmark gene-set comparison using a for loop
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(upreg_1000 $ SYMBOL, 
                          geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000 $ SYMBOL, 
                            geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), 
                                                       up $ p.value, 
                                                       up $ estimate, 'UP')
  hallmarks_results[nrow(hallmarks_results) +1, ] <- c(names(geneid), 
                                                       down $ p.value, 
                                                       down $ estimate, 'DOWN')
}
hallmarks_results <- hallmarks_results %>% 
  mutate(pvalue = as.numeric(pvalue), 
         estimate = as.numeric(estimate))
View(hallmarks_results)

#To store results of fisher test for kegg gene-set
kegg_results <- data.frame(setname = character(), 
                           pvalue = numeric(), estimate = numeric(), 
                           exp = character(), stringsAsFactors = FALSE)

#To store the results for kegg gene-set comparison while using a for loop
for (i in 1:length(KEGG))
{
  geneid <- geneIds(KEGG[i])
  fisher_up <- fishertest(upreg_1000 $ SYMBOL, 
                          geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000 $ SYMBOL, 
                            geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), 
                                             up $ p.value, 
                                             up $ estimate, 'UP')
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid),
                                             down $ p.value, 
                                             down $ estimate, 'DOWN')
}

kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), 
                                        estimate = as.numeric(estimate))
View(kegg_results)

#To store results of fisher test for GO gene-set
go_results <- data.frame(setname = character(), 
                         pvalue = numeric(), estimate = numeric(), 
                         exp = character(), stringsAsFactors = FALSE)

#To store the results for GO gene-set comparison while using a for loop
for (i in 1:length(GO))
{
  geneid <- geneIds(GO[i])
  fisher_up <- fishertest(upreg_1000 $ SYMBOL, 
                          geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(downreg_1000 $ SYMBOL, 
                            geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  go_results[nrow(go_results) +1, ] <- c(names(geneid), 
                                         up $ p.value, 
                                         up $ estimate, 'UP')
  go_results[nrow(go_results) +1, ] <- c(names(geneid), 
                                         down $ p.value, 
                                         down $ estimate, 'DOWN')
}

go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), 
                                    estimate = as.numeric(estimate))
View(go_results)

#adjusting the p-value using benjamini hochberg(BH) method and storing the data in separate file
go_results$BH <- p.adjust(go_results$pvalue, method = "BH",
                          n = length(go_results $ pvalue))
write.csv(go_results, "final_go.csv")     

kegg_results$BH <- p.adjust(kegg_results$pvalue, method = "BH", 
                            n = length(kegg_results $ pvalue))
write.csv(kegg_results, "final_kegg.csv")

hallmarks_results$BH <- p.adjust(hallmarks_results$pvalue, method = "BH", 
                                 n = length(hallmarks_results $ pvalue))
write.csv(hallmarks_results, "final_hallmarks.csv")

# Statistically significant and enriched genesets
# Hallmark:
hallmark_sig <- hallmarks_results[hallmarks_results$pvalue < 0.05,]
hallmark_sig_total <- length(hallmark_sig $ gene_set)
cat("Statistically enriched genesets in Hallmark: ", hallmark_sig_total, "\n")

# Kegg:
kegg_sig <- kegg_results[kegg_results$pvalue < 0.05,]
kegg_sig_total <- length(kegg_sig$gene_set)
cat("Statistically enriched genesets in KEGG: ", kegg_sig_total, "\n")

# GO:
go_sig <- go_results[go_results$pvalue < 0.05,]
go_sig_total <- length(go_sig$gene_set)
cat("Statistically enriched genesets in GO: ", go_sig_total, "\n")

# Top 3 genesets for each:
top3_kegg <- slice_min(kegg_results, order_by = pvalue, n=3)
top3_go <- slice_min(go_results, order_by = pvalue, n=3)
top3_hallmark <- slice_min(hallmarks_results, order_by = pvalue, n=3)
top_3 <- rbind(top3_kegg, top3_go, top3_hallmark)
print(top_3)
write.csv(top_3, file="geneset_final_3_results.csv")