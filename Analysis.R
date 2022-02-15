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
  exp_data <- readr::read_delim(filepath, delim = " ")
  return(exp_data)
}

expr_data <- load_expression('example_intensity_data.csv')


filter_data <- function(tibble){
  
  data_copy <- select(tibble, -probe)
  percent <- 0.20 * (ncol(tibble) -1)
  threshold <- log2(15)
  thresh_tibble <- as_tibble(sapply(data_copy, function(x){ifelse(x>threshold, 1, 0)}))
  count_tibble <- tibble(apply(thresh_tibble, 1, sum))
  combined_tibble <- mutate(tibble, threshold = count_tibble)
  filtered_tibble <- dplyr::filter(combined_tibble, threshold > percent) %>% select(-last_col())
  
  return(filtered_tibble)
}

filter_20 <- filter_data(expr_data)


variance_filter <- function(tibble) {

  probes <- select(tibble, probe)
  data_copy <- select(tibble, -probe)
  df = ncol(data_copy) - 1
  
  chi_lower = qchisq((0.01)/2, df) 
  chi_upper = qchisq((1 - 0.01)/2, df, lower.tail = FALSE)
  
  variance <- apply(data_copy, 1, var)
  test_stat <- (df*variance/median(variance))
  
  filtered_data <- mutate(tibble, stats = test_stat) %>% dplyr::filter(test_stat > chi_upper) %>% select(-stats)
  
  return(filtered_data)
}

filter_variance <- variance_filter(filter_20)


cov_filter <- function(tibble) {
  CVA <- function(x){
    sd(x)/mean(x)}
  
  threshold <- 0.186
  data_copy <- select(tibble, -probe)
  cov_data <- tibble(apply(data_copy, 1, CVA))
  
  filtered_data <- mutate(tibble, coef_var = cov_data) %>% dplyr::filter(coef_var > threshold ) %>% select(-last_col())
  
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
set.seed(240)  # Setting seed
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
annotation_data <- select(annotation_matrix, c("cit-coloncancermolecularsubtype","geo_accession"))

data <- as.matrix(Texpr_data[,-1])
data <-t(data)
data_names <- as.matrix(Texpr_data)

  
metadata <- annotation_data %>% 
  dplyr::right_join(Texpr_data,
  by=c("geo_accession" = "probe")) %>%
  rename("condition" = "cit-coloncancermolecularsubtype")

condition_colors <-
  transmute(
    metadata,
    color=if_else(condition == "C3","red","blue")
  )
condition_colors <- t(condition_colors)

heatmap(data, Colv = NA,
        scale="column",
        cexRow=1,
        cexCol = 0.9,
        ColSideColors=condition_colors,
        labCol=paste(data_names[,1],sep=""))


