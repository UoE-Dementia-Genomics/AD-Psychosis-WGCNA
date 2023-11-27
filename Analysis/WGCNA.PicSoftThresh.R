args<-commandArgs(TRUE)

in_dir <- args[1]
block_size <- as.numeric(args[2])
out_pref <- args[3]

print("Loading libraries...")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
##########################################
#Step 1: Pre-processing 
#########################################
setwd(in_dir)

print("Reading Previous step output...")
if(!file.exists(paste0(out_pref,".wgcna1.rdat"))){
  print("Previous step output doesn't exist!")
}else{
  load(paste0(out_pref,".wgcna1.rdat"))
  print("Phenotype file and Expr file are matched?")
  
  flag=identical(rownames(pheno), colnames(betas))
  if(flag){
     print("yes")
     print("Calculating soft powers...")
     powers = c(c(1:10), seq(from = 12, to=20, by=2))
     sft = pickSoftThreshold(t(betas), powerVector = powers, verbose = 5, dataIsExpr = TRUE,blockSize = block_size)
     save(sft,file = paste0(out_pref,".wgcna.softThreshold.rdat"))
     
     sizeGrWindow(9, 5)
     par(mfrow = c(1,2));
     cex1 = 0.9;
     
     print("Saving plots...")
     pdf(file = paste0(out_pref,".wgcna.softThreshold.pdf"))
     plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
          xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
          main = paste("Scale independence"));
     text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
          labels=powers,cex=cex1,col="red");
     abline(h=0.90,col="red")
     plot(sft$fitIndices[,1], sft$fitIndices[,5],
          xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
          main = paste("Mean connectivity"))
     text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
     dev.off()
     
     print("All done!")
     }else{
          print("NO")
          print("Match phenotype and methylation data before running the script")
     }

}
