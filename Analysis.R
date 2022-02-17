#!/usr/bin/Rscript
## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF528
## Part: Analyst
## Project 1

library(tidyverse)


#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filepath) {
  exp_data <- readr::read_delim(filepath, delim = ",")
  exp_data %>% rownames("probe")
  
  return(exp_data)
}

#expr_data <- load_expression('example_intensity_data.csv')
expr_data <- load_expression('project_1.csv')


filter_data <- function(tibble){
  
  data_copy <- dplyr::select(tibble, -probe)
  percent <- 0.20 * (ncol(tibble) -1)
  threshold <- log2(15)
  thresh_tibble <- as_tibble(sapply(data_copy, function(x){ifelse(x>threshold, 1, 0)}))
  count_tibble <- tibble(apply(thresh_tibble, 1, sum))
  combined_tibble <- dplyr::mutate(tibble, threshold = count_tibble)
  filtered_tibble <- dplyr::filter(combined_tibble, threshold > percent) %>% dplyr::select(-last_col())
  
  return(filtered_tibble)
}

filter_20 <- filter_data(expr_data)


variance_filter <- function(tibble) {
  
  probes <- dplyr::select(tibble, probe)
  data_copy <- dplyr::select(tibble, -probe)
  df = ncol(data_copy) - 1
  
  chi_lower = qchisq((0.01)/2, df) 
  chi_upper = qchisq((1 - 0.01)/2, df, lower.tail = FALSE)
  
  variance <- apply(data_copy, 1, var)
  test_stat <- (df*variance/median(variance))
  
  filtered_data <- dplyr::mutate(tibble, stats = test_stat) %>% dplyr::filter(test_stat > chi_upper) %>% dplyr::select(-stats)
  
  return(filtered_data)
}

filter_variance <- variance_filter(filter_20)


cov_filter <- function(tibble) {
  CVA <- function(x){
    sd(x)/mean(x)}
  
  threshold <- 0.186
  data_copy <- dplyr::select(tibble, -probe)
  cov_data <- tibble(apply(data_copy, 1, CVA))
  
  filtered_data <- dplyr::mutate(tibble, coef_var = cov_data) %>% dplyr::filter(coef_var > threshold ) %>% dplyr::select(-last_col())
  
  return(filtered_data)
}

filter_cov <- cov_filter(filter_variance)

## Write out a different file containing the gene expression matrix for genes 
#' passing all three of the filters from 4.1, 4.2, and 4.3.
## For groups with Biologist role only: Write out the expression matrix for 
#' probesets that pass the expression threshold from 4.2 to a file with write.csv.

## Writting CSV files

write.csv(filter_cov, 
          "/Users/kyragriffin/BU/Spring_2022/BF528/project-1-frazzled-1/gene_expression_matrix.csv", 
          row.names = FALSE)

var_filter <- variance_filter(expr_data)
write.csv(var_filter, 
          "/Users/kyragriffin/BU/Spring_2022/BF528/project-1-frazzled-1/probeset_expression_matrix.csv", 
          row.names = FALSE)

print(paste0("The number of genes that pass all of these thresholds: ", 
             length(filter_cov$probe) - 1))


##TRANSPOSING MATRIX TO CLUSTER

transpose_expression <- function(filtered_expr_data) {
  
  exp_data <- as_tibble(cbind(nms = names(filtered_expr_data), t(filtered_expr_data)))
  colnames(exp_data) <- exp_data[1,]
  exp_data <- exp_data[-1, ]
  
  final_data <- exp_data %>% mutate(across(c(2:1483), as.double))
  
  return (final_data)
}

Texpr_data <- transpose_expression(filter_cov)

##CLUSTERING

# Finding distance matrix
distance_mat <- dist(Texpr_data, method = 'euclidean')
distance_mat

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(112)  # Setting seed
Hierar_cl <- hclust(distance_mat, method = "average")
Hierar_cl

# Plotting dendrogram
plot(Hierar_cl)

# Cutting tree by no. of clusters
fit <- cutree(Hierar_cl, k = 2)
fit

clusters <- table(fit) # Number of samples for each cluster
rect.hclust(Hierar_cl, k = 2, border = "green")

print(paste0("The number of samples in cluster 1 and 2: ", clusters[1]," and ", clusters[2]))


## MAKING HEATMAP

#Getting subtype information
annotation_matrix <- readr::read_delim("proj_metadata.csv", delim = ",")
annotation_data <- dplyr::select(annotation_matrix, c("cit-coloncancermolecularsubtype","geo_accession"))

cond <- dplyr::select(annotation_data, "cit-coloncancermolecularsubtype") %>% 
  dplyr::rename(condition = "cit-coloncancermolecularsubtype")

DF <- dplyr::select(Texpr_data, -1)

DF_data <- as.matrix(sapply(DF, as.numeric)) 
DF_data <- t(DF_data)

data_names <- as.matrix(Texpr_data)


metadata <- annotation_data %>% 
  dplyr::right_join(Texpr_data,
                    by=c("geo_accession" = "probe")) %>% 
  dplyr::rename(condition = "cit-coloncancermolecularsubtype")

metadata <- metadata %>% dplyr::mutate(cond = cond, .before = geo_accession)
metadata <- dplyr::select(metadata, -1)
metadata <- metadata %>% dplyr::rename(condition = cond)

condition_colors <-
  transmute(
    metadata,
    color=if_else(condition == "C3","red","blue")
  )
condition_colors <- t(condition_colors)

heatmap(DF_data,Colv = NA,
        scale="column",
        cexRow=0.85,
        cexCol = 0.35,
        ColSideColors=condition_colors,
        labCol=paste(data_names[,1],sep=""))

####### t test for filtered expression matrix ###############################################
transpose_expr <- function(expr_data) {
  
  exp_data <- as_tibble(cbind(nms = names(expr_data), t(expr_data)))
  colnames(exp_data) <- exp_data[1,]
  exp_data <- exp_data[-1, ]
  
  final_data <- exp_data %>% mutate(across(c(2:ncol(exp_data)), as.double))
  
  return (final_data)
}

## identify genes deferentially expressed between the two clusters using a Welch t-test
fit_tib <- as_tibble(fit) %>% dplyr::rename(cluster = value)
x <- fit_tib %>% dplyr::filter(cluster == 1)
y <- fit_tib %>%  dplyr::filter(cluster == 2)

f_data <- Texpr_data %>%mutate(cluster = fit_tib, .after = probe)
data <- f_data

#calculate_t <- function(data) {
  C3_data <- data %>% 
    dplyr::filter(cluster==1) %>% dplyr::select(!cluster)
  
  C3_table <- transpose_expr(C3_data)
  
  print(ncol(C3_table))
  
  C4_data <- data %>% 
    dplyr::filter(cluster==2)%>%
    dplyr::select(!cluster)
  
  C4_table <- transpose_expr(C4_data)
  
  print(ncol(C4_table))
  
  col_num_C3<-ncol(C3_table)
  col_num_C4<-ncol(C4_table)
  t_value<-c()
  p_value<-c()
  
  
  for (i in 1:nrow(C3_table)){
    temp<-t.test(as.numeric(C3_table[i,1:col_num_C3]),as.numeric(C4_table[i,1:col_num_C4]))
    
    t_value[i]<- temp$statistic
    p_value[i]<- temp$p.value
  }
  probe_name<- colnames(data)
  print(probe_name)
  
  test <- tibble(
    probe = probe_name[3:length(probe_name)],
    t = t_value,
    p_val = p_value) %>%
    arrange(p_val) %>%
    mutate(adj_p = p.adjust(p_val,method="fdr",n=length(p_val)))



count = dplyr::filter(test, adj_p < 0.05) %>% nrow()
print(paste0("The number of differentially expressed genes at adjusted 𝑝<0.05 between the clusters : ", count ))

write.csv(test, 
          "/Users/kyragriffin/BU/Spring_2022/BF528/project-1-frazzled-1/differentially_expressed-genes.csv", 
          row.names = FALSE)

################ t test for probeset expression matrix ###########################################
var_filter <- variance_filter(expr_data)
no_filter_data <- transpose_expression(var_filter)

##CLUSTERING

# Finding distance matrix
distance_mat <- dist(no_filter_data, method = 'euclidean')
distance_mat

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(112)  # Setting seed
Hierar_cl <- hclust(distance_mat, method = "average")
Hierar_cl

# Plotting dendrogram
plot(Hierar_cl)

# Cutting tree by no. of clusters
fit <- cutree(Hierar_cl, k = 2)
fit

fit_tib <- as_tibble(fit) %>% dplyr::rename(cluster = value)
x <- fit_tib %>%  dplyr::filter(cluster == 1)
y <- fit_tib %>%  dplyr::filter(cluster == 2)

no_filter_data <- no_filter_data %>% dplyr::mutate(cluster = fit_tib, .after = probe)


  probe_C3_data <- no_filter_data %>% 
    dplyr::filter(cluster==1) %>% dplyr::select(!cluster)
  
  probe_C3_table <- transpose_expr(probe_C3_data)
  
  print(ncol(probe_C3_table))
  
  probe_C4_data <- no_filter_data %>%
    dplyr::filter(cluster==2)%>%
    dplyr::select(!cluster)
  
  probe_C4_table <- transpose_expr(probe_C4_data)
  
  print(ncol(probe_C4_table))
  
  col_num_C3<-ncol(probe_C3_table)
  col_num_C4<-ncol(probe_C4_table)
  t_value<-c()
  p_value<-c()
  
  
  for (i in 1:nrow(probe_C3_table)){
    temp<-t.test(as.numeric(probe_C3_table[i,1:col_num_C3]),as.numeric(probe_C4_table[i,1:col_num_C4]))
    
    t_value[i]<- temp$statistic
    p_value[i]<- temp$p.value
  }
  probe_name<- colnames(no_filter_data)
  print(probe_name)
  
  probe_test <- tibble(
    probe = probe_name[3:length(probe_name)],
    t = t_value,
    p_val = p_value) %>%
    arrange(p_val) %>%
    mutate(adj_p = p.adjust(p_val,method="fdr",n=length(p_val)))%>%
    filter(adj_p < 0.05)

  

count = dplyr::filter(probe_test, adj_p < 0.05) %>% nrow()
print(paste0("The number of differentially expressed genes at adjusted 𝑝<0.05 between the clusters : ", count ))

write.csv(probe_test, 
          "/Users/kyragriffin/BU/Spring_2022/BF528/project-1-frazzled-1/probeset_differentially_expressed-genes.csv", 
          row.names = FALSE)  
  
  
  
  
  
  
  


