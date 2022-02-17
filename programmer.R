#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

#BiocManager::install("affy"")
#BiocManager::install("affyPLM")
#BiocManager::install("sva")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("hgu133plus2.db")
install.packages("ggfortify")

#module load R
#R --vanilla

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggfortify)


dataset1 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/CEL_files/')
dataset2 <- ReadAffy(celfile.path = '/projectnb/bf528/users/frazzled/project_1/samples/')


files2<-merge(dataset1,dataset2) 
files2

rma_data <- affy::rma(files2)
fit_data <- fitPLM(files2,normalize=TRUE,background=TRUE)
RLE_data <-RLE(fit_data,type="stats")
NUSE_data <- NUSE(fit_data,type="stats")

RLE_hist <- hist(RLE_data[1,])
NUSE_hist <- hist(NUSE_data[1,])

metadata <- readr::read_csv('/project/bf528/project_1/doc/proj_metadata.csv')
batch <- metadata$normalizationcombatbatch
mod <- model.matrix(~as.factor(normalizationcombatmod), data = metadata)

edata = (exprs(rma_data))
combat_data = ComBat(edata,batch,mod)

write.csv(combat_data, '/projectnb/bf528/users/frazzled/project_1.csv')

transposed_data <- t(edata)
scaled_data <- scale(transposed_data)
data_for_pca <- t(scaled_data)

pca_results <- prcomp(data_for_pca, center=FALSE, scale=FALSE)

pca <- pca_results$rotation

res <- summary(pca_results)
res2 <- summary(pca)


PC1_imp = cat("PC1 percent variability = ",res$importance[2,1])
PC2_imp = cat("PC2 percent variability = ",res$importance[2,2])


plot(pca_results$rotation[,1:2], xlab = PC1_imp, ylab = PC2_imp)

