args <- commandArgs(TRUE)

beta_file <- args[1] 
pheno_file <- args[2]
trait <- args[3]
variables_fact <- args[4]
variables_num <- args[5]
mad_thr <- as.numeric(args[6])
out_pref <- args[7]

sink(paste0(out_pref,".wgcna.Preprocessing.log.txt"))

print("Input arguments:")
print(paste("    Beta value file:",beta_file))
print(paste("    Phenotype file:",pheno_file))
print(paste("    Trait variable:",trait))
print(paste("    Factor variables:",variables_fact))
print(paste("    Numeric variables:",variables_num))
print(paste("    MAD Threshold:",mad_thr))
print(paste("    Output prefix:",out_pref))
cat("\n")


print("Loading libraries...")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dendextend))

##########################################
#Step 1: Pre-processing 
#########################################

print("Reading phenotypes and beta values...")
betas <- readRDS(file=beta_file)
pheno = read.csv(pheno_file,stringsAsFactors = F,row.names = 1)

print("Phenotype file and betas file are matched?")
print(identical(rownames(pheno),colnames(betas)))

print(paste("Keeping top",(mad_thr*100),"% of CpGs with high Median absolute deviation"))
## calculate gene Median Absolute Deviation across all samples
gene_mad <- apply(betas,1,mad)

gene_mad <- sort(gene_mad,decreasing = T)
gene_mad1 <- gene_mad[1:round(length(gene_mad)*mad_thr,digits = 0)]

print(paste("Total number of CpGs:",length(gene_mad)))
print(paste("Number of CpGs after filtering:",length(gene_mad1)))
## removing genes with low MAD
betas <- betas[rownames(betas) %in% rownames(as.data.frame(gene_mad1)),]

## checking genes and samples with too many missing value
passed_gene_samples <- goodSamplesGenes(betas,verbose = 3)
print("goodSamplesGenes function: All samples and genes passed?")
passed_gene_samples$allOK
# [1] True

## finding outlier samples
print("Hclust...")

distance <- dist(t(betas),method = "euclidean")
samples_tree <- hclust(distance,method = "average")

variables_num <- str_split(variables_num,pattern = ',',simplify = T)[1,]
variables_fact <- str_split(variables_fact,pattern = ',',simplify = T)[1,]

pheno2 <- as.numeric(as.factor(pheno[,trait]))
for (i in 1:length(variables_fact)) {
  pheno2 <- cbind.data.frame(pheno2,as.numeric(as.factor(pheno[,variables_fact[i]])))
}
for (i in 1:length(variables_num)) {
  pheno2 <- cbind.data.frame(pheno2,as.numeric(pheno[,variables_num[i]]))
}
names(pheno2) <- c(trait,variables_fact,variables_num)

traitColors = numbers2colors(pheno2, signed = FALSE)

pdf(file = paste0(out_pref,".hclust.pdf"),width = 20,height = 20)
par(cex = 0.6)
par(mar = c(0,4,2,0)) 
plot(samples_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
## Plot the sample dendrogram and the colors underneath. 
plotDendroAndColors(samples_tree, traitColors, groupLabels = names(pheno2), main = "Sample dendrogram and trait heatmap")
dev.off()

print("Calculating PCs...")
pca <- prcomp(t(betas),center = T,scale. = T)
pca1 <- pca$x[,1:10]
pca1 <- data.frame(pca1)

cor_ <- matrix(data = NA,nrow =10,ncol = length(c(trait,variables_num,variables_fact)) )
colnames(cor_) <- c(trait,variables_num,variables_fact)
rownames(cor_) <- colnames(pca1)[1:10]

cor_pval <- cor_
print("Calculatin correlations...")

for (i in 1:10) {
  for (j in 1:length(c(trait,variables_num,variables_fact))) {
    var_ <- c(trait,variables_num,variables_fact)[j]
    
    if(var_ %in% c(trait,variables_fact)){
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(as.factor(pheno2[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(pheno2[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)
print("Saving Plot...")
pdf(file = paste0(out_pref,".PCA.pdf"))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = rownames(cor_),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))

dev.off()

print("Saving data...")
save(betas, pheno, file = paste0(out_pref,".wgcna1.rdat"))

print("All done!")

sink()
