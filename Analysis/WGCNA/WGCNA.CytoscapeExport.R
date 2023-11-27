args<-commandArgs(TRUE)

net_file <- args[1]
module_name <- args[2]
weighted_ <- ifelse(args[3]=="T",T,F)
threshold_ <- as.numeric(args[4])
out_pref <- args[5]

print("Input arguments:")
print(paste("    Network file:",net_file))
print(paste("    Module name:",module_name))
print(paste("    Is weighted?",weighted_))
print(paste("    Adjacency threshold:",threshold_))
print(paste("    Output prefix:",out_pref))
cat("\n")
print("Loading libraries...")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

#########################################

print("Openning network file...")
net = readRDS(net_file)

print("Finding Module nodes and blocks...")
module.nodes.index = which(net$colors==module_name)
module.blocks = net$blocks[module.nodes.index]

print("Module are in the following blocks:")
print(unique(module.blocks))

for (b in unique(module.blocks)) {
  
  print(paste("Loading block",b,"..."))
  load(paste0(dirname(net_file),"/",net$TOMFiles[b]))
  
  print(paste("Convert block",b,"To matrix ..."))
  TOM.mat = as.matrix(TOM)
  module.nodes.TOM.index = net$blockGenes[[b]] %in% module.nodes.index
  module.TOM.mat=TOM.mat[module.nodes.TOM.index,module.nodes.TOM.index]
  
  print(paste("Export Cytoscape data for block",b,"..."))
  cytoscape = exportNetworkToCytoscape(adjMat = module.TOM.mat,weighted = weighted_ , threshold = threshold_,
                                       edgeFile = paste0(out_pref,".",module_name,".Block.",b,".Thresh.",threshold_,".egdeData.txt"),
                                       nodeFile = paste0(out_pref,".",module_name,".Block.",b,".Thresh.",threshold_,".nodeData.txt"),
                                       nodeNames = names(module.nodes.index[module.blocks==b]))
}
