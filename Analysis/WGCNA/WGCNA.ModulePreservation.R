library(WGCNA)
allowWGCNAThreads()
library(stringr)

args <- commandArgs(T)

exp_file_1 <- args[1]
exp_file_2 <- args[2]
net_file_1 <- args[3]
modules <- args[4]
max_size <- as.numeric(args[5])
max_gold_size <- as.numeric(args[6])
n_permut <- as.numeric(args[7])
out_pref <- args[8]

print(paste0("Maximum module size set to : ",max_size))

exp1 <- readRDS(exp_file_1)
exp2 <- readRDS(exp_file_2)
net1 <- readRDS(net_file_1)
modules <- str_split(modules,pattern = ',',simplify = T)[1,]

selected.colors = net1$colors[net1$colors %in% modules]
selected.exp1 = exp1[names(selected.colors),]

multiExpr  = list(A1=list(data=t(selected.exp1)),A2=list(data=t(exp2))) 
multiColor = list(A1 = selected.colors) 

mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,nPermutations = n_permut, networkType="unsigned",maxModuleSize = max_size,parallelCalculation = F,maxGoldModuleSize = max_gold_size,savePermutedStatistics = F )
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2  
save(mp,stats,file = paste0(out_pref,"_ModulePreservation.Rdata"))
print(stats[order(-stats[,2]),c(1:2)])
