#install packages
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")
install.packages("ggplot2")
install.packages("ggpubr")
install_github("vqv/ggbiplot")

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(ggfortify)

#Read CEL files using ReadAffy to create an affybatch object
dataset1 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/CEL_files/')
dataset2 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/')

#Merge the datasets
files2<-merge(dataset1,dataset2) 

#RMA normalized expression values for the affybatch object
rma_data <- affy::rma(files2)

#Converting an AffyBatch object into an PLMset using fitPLM
fit_data <- fitPLM(files2,normalize=TRUE,background=TRUE)

#Compute RLE (relative log expression) of the microarray samples
RLE_data <-RLE(fit_data,type="stats")

#Compute NUSE (noramalized unscaled standard error)
NUSE_data <- NUSE(fit_data,type="stats")

#Plotting the medians of RLE and NUSE in histograms
RLE_hist <- hist(RLE_data[1,], xlab = "Median RLE", col = "lightgreen", main = NULL)
NUSE_hist <- hist(NUSE_data[1,], xlab = "Median NUSE", col = "lightblue", main =NULL)


#Path to the annotation file on SCC
metadata <- readr::read_csv('/project/bf528/project_1/doc/proj_metadata.csv')

# Batch effects variable
batch <- metadata$normalizationcombatbatch

#Creating a model matrix containing covariates
mod <- model.matrix(~as.factor(normalizationcombatmod), data = metadata)

#Transforming normalized data into a matrix
edata = (exprs(rma_data))

#Correcting batch effects while preserving features of interest
combat_data = ComBat(edata,batch,mod)

#Writing the ComBat file
write.csv(combat_data, '/projectnb/bf528/users/frazzled/samples/project_1.csv', ifelse(append, "a", "w"))

# Transforming, scaling and retransposing the data 
transposed_data <- t(edata)
scaled_data <- scale(transposed_data)
data_for_pca <- t(scaled_data)

#perform a principal components analysis on the given matrix
pca_results <- prcomp(data_for_pca, center=FALSE, scale=FALSE)

#Creating a data frame of PC1 and PC2 components
pc12 <- as.data.frame(pca_results$rotation[,1:2])

#Plot PC1 and PC2 with variance
ggplot(data = pc12, aes(PC1,PC2)) + geom_point(color='red')+
  theme_light()+ xlab('PC1 (14.5%)')+ylab('PC2 (9.54%)')+ggtitle('PCA of Gene Expression for 134 samples')

#Identifying the outliers with a box plot
boxplot(pc12$PC1, pc12$PC2, main = 'Boxplots for PC1 and PC2', names = c("PC1","PC2"), col = "lightblue", border = "purple", notch = TRUE)

#Creating a data frame of the pca after using $rotation method
pca <- as.data.frame(pca_results$rotation)

#Trying to count the samples that are 3 standard deviations away from the mean
pc_out <- which(pca$PC1 > mean(pca$PC1) + 3*sd(pca$PC1) | pca$PC1 < mean(pca$PC1) - 3*sd(pca$PC1) |
                  pca$PC2 > mean(pca$PC2) + 3*sd(pca$PC2) | pca$PC2 < mean(pca$PC2) - 3*sd(pca$PC2))
(pc_out)

