Module.Membership.Plot <- function(net.colors, expr.mat, trait, modules, size.threshold,soft.power,plot.title=NA,selected.genes=NA){
  suppressMessages(library(WGCNA))
  suppressMessages(library(ggplot2))
  
  if(!identical(names(net.colors),rownames(expr.mat))){
    index <- intersect(names(net.colors),rownames(expr.mat))
    net.colors <- net.colors[index]
    expr.mat <- expr.mat[index, ]
  }
  print("Calculating Modules Eigengenes")
  MEs = moduleEigengenes(expr = t(expr.mat) , colors = net.colors,softPower = soft.power)$eigengenes
  MEs <- orderMEs(MEs)
  names(MEs) = substring(names(MEs), 3)
  MEs <- MEs[,modules]
  nSamples = ncol(expr.mat)
  
  print("Calculating Modules Memberships and P values...")
  geneModuleMembership = as.data.frame(cor(t(expr.mat), MEs, use = "p"))
  names(geneModuleMembership) <- modules
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  MMPvalue.adj = as.data.frame(apply(MMPvalue,2,function(x){p.adjust(x,method = "bonferroni")}))
  
  print("Calculating Gene Significant and P values...")
  geneTraitSignificance = as.data.frame(cor(t(expr.mat), as.numeric(as.factor(trait[,1])), method = "spearman"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  GSPvalue.adj = apply(GSPvalue,2,function(x){p.adjust(x,method = "bonferroni")})
  
  names(geneTraitSignificance) = names(trait)
  names(GSPvalue) = names(trait)
  
  Plots <- vector(mode = "list",length = length(modules))
  Corrs <- data.frame(Module=vector(mode = "character",length = length(modules)),
                      Size=vector(mode = "character",length = length(modules)),
                      Cor=vector(mode = "character",length = length(modules)),
                      Pvalue=vector(mode = "character",length = length(modules)))
  names(Plots) <- modules
  for (i in 1:length(modules)) {
    moduleGenes = (net.colors==modules[i])
    
    if(length(moduleGenes[moduleGenes])>size.threshold){
      mm = abs(geneModuleMembership[moduleGenes , modules[i]])
      mm = order(mm,decreasing = T)
      moduleGenes = rownames(geneModuleMembership)[1:round(length(mm)*0.2,digits = 0)]
    }
    plot.data <- cbind.data.frame(ID = rownames(geneModuleMembership)[moduleGenes], MM=abs(geneModuleMembership[moduleGenes , modules[i]]),GS=abs(geneTraitSignificance[moduleGenes , 1]))
    c <- cor.test(plot.data$MM,plot.data$GS,method = "pearson")
    print(paste("Generating Plot for",modules[i],"module..."))
    p <- ggplot(data = plot.data,aes(x=MM,y=GS))+
      geom_point(color=modules[i])+
      geom_smooth(method = "lm") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))+
      xlab("Module Membership") +
      ylab("Probe Significance") 
    
    if(is.na(plot.title)){
      p <- p+ggtitle(paste0(modules[i],"\n","(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }else{
      p <- p+ggtitle(paste0(plot.title,"\n","(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }
    if(!is.na(selected.genes[1])){
      p <- p + geom_point(data=plot.data[plot.data$ID %in% selected.genes,],
                          pch=21, fill="blue", size=4, colour="blue", stroke=1)
    }
    Corrs$Module[i] <- modules[i]
    Corrs$Size[i] <- length(net.colors[net.colors==modules[i]])
    Corrs$Cor[i] <- c$estimate
    Corrs$Pvalue[i] <- c$p.value
    Plots[[i]] <- p
  }
  return(list(Plots = Plots , Cor.test = Corrs))
}
