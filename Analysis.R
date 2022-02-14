#!/usr/bin/Rscript
## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF528
## Part: Analyst
## Project 1

library(tidyverse)

#' #### Loading and processing data ####
#' #' Load Expression Data
#' #'
#' #' @param filepath A text string of the full filepath to the file to load.
#' #'
#' #' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' #' 
#' #' @details Note that not all CSVs are created equal, and there are often cases where 
#' #' the data will not load in correctly on the first try. You may want to write this functon to 
#' #' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' #' tibble correctly.
#' #'
#' #' @examples 
#' #' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
#' load_expression <- function(filepath) {
#'   
#'   exp_data <- t(readr::read_delim(filepath, delim = " ", col_names = TRUE))
#'   colnames(exp_data) <- exp_data[1,]
#'   exp_data <- exp_data[-1, ] 
#'   
#'   exp_data <- as_tibble(exp_data, rownames = "subject_id")
#'   final_data <- exp_data %>% mutate(across(c(2:54676), as.double))
#'   
#'   return (final_data)
#' }

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
