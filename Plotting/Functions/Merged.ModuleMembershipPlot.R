Merged.Module.Membership.Plot <- function(net.colors, expr.mat1, expr.mat2, trait1, trait2, modules, 
                                   soft.power, plot.title,legend.title,legend.value){
  suppressMessages(library(WGCNA))
  suppressMessages(library(ggplot2))
  
  index <- intersect(names(net.colors),rownames(expr.mat1))
  net.colors1 <- net.colors[index]
  expr.mat1 <- expr.mat1[index, ]
  
  index <- intersect(names(net.colors),rownames(expr.mat2))
  net.colors2 <- net.colors[index]
  expr.mat2 <- expr.mat2[index, ]
  
  ################################################################################
  print("Calculating Modules Eigengenes based on expression matrix 1...")
  MEs1 = moduleEigengenes(expr = t(expr.mat1) , colors = net.colors1,softPower = soft.power)$eigengenes
  print("Calculating Modules Eigengenes based on expression matrix 2...")
  MEs2 = moduleEigengenes(expr = t(expr.mat2) , colors = net.colors2,softPower = soft.power)$eigengenes
    
  MEs1 <- orderMEs(MEs1)
  MEs2 <- orderMEs(MEs2)
  names(MEs1) <- substring(names(MEs1), 3)
  names(MEs2) <- substring(names(MEs2), 3)
  MEs1 <- MEs1[,modules]
  MEs2 <- MEs2[,modules]
 
  nSamples1 = ncol(expr.mat1)
  nSamples2 = ncol(expr.mat2)
  #############################################################################
  print("Calculating Modules Memberships and P values based on expression matrix 1...")
  geneModuleMembership1 = as.data.frame(cor(t(expr.mat1), MEs1, use = "p"))
  names(geneModuleMembership1) <- modules
  MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples1));
  MMPvalue1.adj = as.data.frame(apply(MMPvalue1,2,function(x){p.adjust(x,method = "bonferroni")}))
  
  ##############################################################################
  print("Calculating Modules Memberships and P values based on expression matrix 2...")
  geneModuleMembership2 = as.data.frame(cor(t(expr.mat2), MEs2, use = "p"))
  names(geneModuleMembership2) <- modules
  MMPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership2), nSamples2));
  MMPvalue2.adj = as.data.frame(apply(MMPvalue2,2,function(x){p.adjust(x,method = "bonferroni")}))
  
  ##################################################################################
  print("Calculating Gene Significant and P values based on expression matrix 1...")
  geneTraitSignificance1 = as.data.frame(cor(t(expr.mat1), as.numeric(as.factor(trait1[,1])), method = "spearman"))
  GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples1))
  GSPvalue1.adj = apply(GSPvalue1,2,function(x){p.adjust(x,method = "bonferroni")})
  
  names(geneTraitSignificance1) = names(trait1)
  names(GSPvalue1) = names(trait1)
  ###################################################################################
  print("Calculating Gene Significant and P values based on expression matrix 2...")
  geneTraitSignificance2 = as.data.frame(cor(t(expr.mat2), as.numeric(as.factor(trait2[,1])), method = "spearman"))
  GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples2))
  GSPvalue2.adj = apply(GSPvalue2,2,function(x){p.adjust(x,method = "bonferroni")})
  
  names(geneTraitSignificance2) = names(trait2)
  names(GSPvalue2) = names(trait2)
  
  ################################################################################
  result <- vector(mode = "list",length = length(modules))
  names(result) <- modules
  for (i in 1:length(modules)) {
    moduleGenes1 = (net.colors1==modules[i])
    moduleGenes2 = (net.colors2==modules[i])
    
    plot.data1 <- cbind.data.frame(x=abs(geneModuleMembership1[moduleGenes1 , modules[i]]),y=abs(geneTraitSignificance1[moduleGenes1 , 1]))
    c1 <- cor.test(plot.data1$x,plot.data1$y,method = "pearson")
    
    plot.data2 <- cbind.data.frame(x=abs(geneModuleMembership2[moduleGenes2 , modules[i]]),y=abs(geneTraitSignificance2[moduleGenes2 , 1]))
    c2 <- cor.test(plot.data2$x,plot.data2$y,method = "pearson")
    
    print(paste("Generating Plot for",modules[i],"module..."))
    p <- ggplot(data =plot.data1 , aes(y = y, x = x)) + ggplot2::geom_point(color="steelblue") +
          geom_smooth( data = plot.data1, aes(y = y, x = x, colour=legend.value[1]), method = "lm",se = T) +
          geom_point(data = plot.data2,aes(y = y, x = x), color = "darkred") +
          geom_smooth(data = plot.data2, aes(y = y, x = x, colour=legend.value[2]),method = "lm",se = T) +
          scale_colour_manual(name=legend.title, breaks =legend.value ,values=c("steelblue","darkred")) + 
          xlab(paste("Module Membership in", modules[i], "module")) +
          ylab(paste("Probe significance for",names(trait1))) + 
          geom_text(x = 0, y = max(c(plot.data1$y,plot.data2$y))+0.01, label = paste(sep = "", "corr = ", round(c1$estimate, 2), ", p = ", 
                    formatC(c1$p.value, format = "e", digits = 2)),color = "steelblue", hjust = 0, vjust = 1) + 
          geom_text(x = 0, y = max(c(plot.data1$y,plot.data2$y)), label = paste(sep = "", "corr = ", round(c2$estimate, 2), ", p = ", 
                   formatC(c2$p.value, format = "e", digits = 2)),color = "darkred", hjust = 0, vjust = 1)+
          theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),legend.key = element_rect(fill = "white", colour = "white"))
      
    
    result[[i]] <- p
    
  }
  return(result)
}
