args<-commandArgs(TRUE)

in_pref <- args[1]
soft_pow <- as.numeric(args[2])
block_size <- as.numeric(args[3])
module_size <- as.numeric(args[4])
out_pref <- args[5]

sink(paste0(out_pref,".wgcna.BlocwiseNet.log.txt"))
print("Input arguments:")
print(paste("    Input prefix:",in_pref))
print(paste("    Soft Power threshold:",soft_pow))
print(paste("    Maximum block size:",block_size))
print(paste("    Minimum module size:",module_size))
print(paste("    Output prefix:",out_pref))
cat("\n")
print("Loading libraries...")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

#########################################

print("Reading Previous step output...")
if(!file.exists(paste0(in_pref,".wgcna1.rdat"))){
  print("Previous step output doesn't exist!")
}else{
  load(paste0(in_pref,".wgcna1.rdat"))
  
  print("Generating network...")
  net = blockwiseModules(t(betas), maxBlockSize = block_size, power = soft_pow, TOMType = "unsigned", minModuleSize = module_size, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = F, saveTOMs= T, verbose = 3, nThreads=16)
  print("Saving the network...")
  save(net,file = paste0(out_pref,".wgcna.network.rdat"))
  
  print("Generating plots...")
  sizeGrWindow(6,6)
  pdf(paste0(out_pref,".wgcna.netDendrogram.pdf"))
  for (i in 1:length(net$dendrograms)) {
    plotDendroAndColors(net$dendrograms[[i]], net$blockGenes[[i]],
                        "Module colors", main = "Gene dendrogram and module colors in block 1", 
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
  dev.off()
  print("All done!")
}

sink()
