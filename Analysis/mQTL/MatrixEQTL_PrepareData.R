args<-commandArgs(trailingOnly = TRUE)

suppressMessages(library(stringr))
suppressMessages(library(data.table))

InDir<-args[1] 
OutDir <- args[2] 
covar_fact <- args[3] #"Sex,TissueType,Plate,BraakStage"
covar_num <- args[4] #"Age,CellProportion"
FilePrefix <- args[5] 

setwd(InDir)

samples <- read.table(file=paste0(FilePrefix,".fam"),sep=" ",header=F,stringsAsFactors=F)
rownames(samples) <- samples$V2

print("Reading genotype eigenvectors...")
eigenvec <- read.table(file=paste0(FilePrefix,".eigenvec"),header=F,row.names = 2)[,c(-1)]

print("Reading expression data...")
exp_all <- as.data.frame(readRDS(paste0(FilePrefix,"_expression.rds")))
exp_sample <- read.csv(paste0(FilePrefix,"_expression.csv"),row.names=1,stringsAsFactors = F) 

if(!all(sapply(list(colnames(exp_all), rownames(exp_sample),rownames(eigenvec)), FUN = identical, rownames(samples)))){
  shared_names <- Reduce(intersect, list(colnames(exp_all),rownames(exp_sample),rownames(samples),rownames(eigenvec)))
  exp_all <- exp_all[,colnames(exp_all) %in% shared_names]
  exp_sample <- exp_sample[rownames(exp_sample) %in% shared_names,]
  samples <- samples[rownames(samples) %in% shared_names,]
  
  exp_all <- exp_all[shared_names]
  exp_sample <- exp_sample[match(shared_names,rownames(exp_sample)),]
  samples <- samples[match(shared_names,rownames(samples)),]
  eigenvec <- eigenvec[match(shared_names,rownames(eigenvec)),]
}

snps_all=vector(mode = "list",length = 22)
snps_loc=vector(mode = "list",length = 22)
for (i in 1:22) {
  print(paste("Reading Genotype data chr",i,"..."))
  if(!file.exists(paste0(FilePrefix,"_chr",i,"_genotype.txt"))){
  
    snps_all[[i]]<-fread(paste0(FilePrefix,"_chr",i,".raw"), data.table = FALSE,nThread=16)
    
    snps_all[[i]] <- t(snps_all[[i]])
    snps_all[[i]] <- cbind.data.frame(rownames(snps_all[[i]]),snps_all[[i]])
    names(snps_all[[i]]) <- c("snpid",snps_all[[i]][2,c(-1)])
    snps_all[[i]] <- snps_all[[i]][c(-1:-6),]
    write.table(snps_all[[i]],file=paste0(FilePrefix,"_chr",i,"_genotype.txt"),sep="\t",row.names = F,col.names = T,quote = F)
  } else {
    snps_all[[i]]<-fread(paste0(FilePrefix,"_chr",i,"_genotype.txt"), data.table = FALSE,nThread=16)
  }
  if(!identical(colnames(exp_all),colnames(snps_all[[i]])[-1])){
    index <- match(colnames(exp_all),colnames(snps_all[[i]])[-1])
    snps_all[[i]] <- snps_all[[i]][,c(1,index+1)]
  }
  snps_loc[[i]] <- str_split(snps_all[[i]]$snpid,pattern=':',simplify=T)[,c(1,2)]
  snps_loc[[i]] <- cbind.data.frame(snps_all[[i]]$snpid,snps_loc[[i]])
  names(snps_loc[[i]]) <- c("snpid","chr","pos")
} 

print("Reading Gene location data...")
gene_loc_all <- read.csv(file=paste0(FilePrefix,"_GeneLocation.csv"),stringsAsFactors = F,row.names = 1)
names(gene_loc_all) <- tolower(names(gene_loc_all))
gene_loc_all <- cbind.data.frame(rownames(gene_loc_all),gene_loc_all$chr,gene_loc_all$start,gene_loc_all$end)
names(gene_loc_all) <- c("geneid","chr","left","right")

if(!identical(row.names(exp_all) , gene_loc_all$geneid)){
  shared_names <- Reduce(intersect, list(row.names(exp_all),gene_loc_all$geneid))
  gene_loc_all <- gene_loc_all[gene_loc_all$geneid %in% shared_names,]
  exp_all <- exp_all[row.names(exp_all) %in% shared_names,]
}

exp_all <- cbind.data.frame(gene_loc_all$geneid,exp_all)
names(exp_all)[1] <- "geneid"

print(paste("Are the gene ids matched in gene expression and location data?"))
print(ifelse(identical(exp_all$geneid,gene_loc_all$geneid),"Yes","NO"))

print("Saving Expression data...")
write.table(exp_all,file = paste0(OutDir,"/",FilePrefix,"_exp.txt"),quote = F,col.names = T,row.names = F,sep = '\t')

print("Saving gene location data...")
write.table(gene_loc_all,file = paste0(OutDir,"/",FilePrefix,"_gene_loc.txt"),quote = F,col.names = T,row.names = F,sep = '\t')

for (i in 1:22) {
  
  print(paste("Are the SNP ids matched in in chr",i," data?"))
  print(ifelse(identical(snps_all[[i]]$snpid,snps_loc[[i]]$snpid),"Yes","NO"))
  
  print(paste("Are the sample ids matched in expression data and chr",i," SNP data?"))
  print(ifelse(identical(colnames(exp_all)[-1], colnames(snps_all[[i]])[-1]),"Yes","NO"))
  
  print(paste("Saving Chr",i,"..."))
  write.table(snps_all[[i]],file = paste0(OutDir,"/",FilePrefix,"_snps","_chr",i,".txt"),quote = F,col.names = T,row.names = F,sep = '\t')
  write.table(snps_loc[[i]],file = paste0(OutDir,"/",FilePrefix,"_snps_loc","_chr",i,".txt"),quote = F,col.names = T,row.names = F,sep = '\t')

}

print("Generating covariates file...")
covar_fact = str_split(covar_fact,pattern = ',',simplify = T)[1,]
covar_num = str_split(covar_num,pattern = ',',simplify = T)[1,]
covariates <- cbind.data.frame(exp_sample[,covar_fact],exp_sample[,covar_num],eigenvec$V3,eigenvec$V4,eigenvec$V5)
names(covariates) <- c(covar_fact,covar_num,"PC1","PC2","PC3")
rownames(covariates) <- rownames(exp_sample)
for (c in covar_fact) {
  covariates[,c] = as.numeric(as.factor(covariates[,c]))
}
covariates <- t(covariates)
covariates <- cbind.data.frame(rownames(covariates),covariates)
names(covariates)[1] <- "id"

print("Sample ids matched in covariate file with the others?")
print(ifelse(all(sapply(list(colnames(exp_all)[-1], rownames(exp_sample),rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"))
#True

print("Saving covaraites file...")
write.table(covariates,file = paste0(OutDir,"/",FilePrefix,"_covariates.txt"),quote = F,col.names = T,row.names = F,sep = '\t')

print("All done!")


