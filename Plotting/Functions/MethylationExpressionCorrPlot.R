methylation.expression.corr <- function(data.methylation,data.expression,correlation.method="pearson",cpg.gene.list,phenotype,No.Pairs.in.Plot=10, trait,trait.labels,return.csv){
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(patchwork)
  if(!identical(rownames(phenotype),colnames(data.expression)) | !identical(rownames(phenotype),colnames(data.methylation))){
    index <- intersect(rownames(phenotype),colnames(data.methylation))
    index <- intersect(index , colnames(data.expression))
    phenotype <- phenotype[index,]
    data.methylation <- data.methylation[,index]
    data.expression <- data.expression[,index]
  }
  cpg.gene.list$Correlation <- NA
  cpg.gene.list$Pvalue <- NA
  
  index <- intersect(cpg.gene.list$cpg , rownames(data.methylation))
  data.methylation <- data.methylation [index , ]
  cpg.gene.list <- cpg.gene.list[cpg.gene.list$cpg %in% index , ]
  
  index <- intersect(cpg.gene.list$gene.ensembl , rownames(data.expression))
  data.expression <- data.expression[unique(cpg.gene.list$gene.ensembl),]
  cpg.gene.list <- cpg.gene.list[cpg.gene.list$gene.ensembl %in% index , ]
  
  for(i in 1:nrow(cpg.gene.list)){
    temp <- cor.test(x = as.numeric(data.methylation[cpg.gene.list$cpg[i],]),y=as.numeric(data.expression[cpg.gene.list$gene.ensembl[i], ]),method =correlation.method)
    cpg.gene.list$Correlation[i] <- temp$estimate
    cpg.gene.list$Pvalue[i] <- temp$p.value
  }
  cpg.gene.list$Padj <- p.adjust(cpg.gene.list$Pvalue,method = "BH")
  cpg.gene.list <- cpg.gene.list[order(cpg.gene.list$Pvalue,decreasing = F),]
  
  result.plot <- vector(mode = "list",length = 10)
  names(result.plot) <- paste0(cpg.gene.list$cpg[1:10],".",cpg.gene.list$gene.symbol[1:10])
  for (i in 1:No.Pairs.in.Plot) {
    data = cbind.data.frame(as.numeric(data.methylation[cpg.gene.list$cpg[i],]),as.numeric(data.expression[cpg.gene.list$gene.ensembl[i],]))
    names(data) = c("Methylation","Expression")
    rownames(data) = colnames(data.expression)
    data$trait = as.factor(phenotype[,trait])
  
  p1 <- ggplot(data = data, aes(x = Methylation, y = Expression)) + geom_point() + 
      geom_smooth(method = "lm", se = T) + theme_bw()+
      ggtitle("")+xlab("Methylation \n(normalized beta value)")+
      ylab("Expression \n(logCPM)")+
     theme(axis.title = element_text(size = 8))
  
   p2 <- ggplot(data = data, aes(x=trait, y=Methylation,fill=trait)) +
     geom_boxplot() +
     guides(fill=guide_legend(title=trait))+
     scale_fill_manual(values = c("#01BEC3","#F8756D"),labels=trait.labels)+
     geom_jitter(color="black", size=0.8, alpha=0.9) +
     ggtitle("") +
     xlab(trait)+ylab("Methylation \n(normalized beta value)") +
     scale_x_discrete(labels=trait.labels)+
     theme_bw()+theme(axis.title = element_text(size = 8))
  
   p3 <- ggplot(data = data, aes(x=trait, y=Expression,fill=trait)) +
     geom_boxplot() +
     guides(fill=guide_legend(title=trait))+
     scale_fill_manual(values = c("#01BEC3","#F8756D"),labels=trait.labels)+
     geom_jitter(color="black", size=0.8, alpha=0.9) +
     ggtitle("") +
     xlab(trait)+ylab("Expression \n(logCPM)") +
     scale_x_discrete(labels=trait.labels)+
     theme_bw()+theme(axis.title = element_text(size = 8))
   
   p4 <- p2+p3+p1+plot_annotation(title = paste0(cpg.gene.list$cpg[i],"-",cpg.gene.list$gene.symbol[i],
                                                 " (Correlation=",round(cpg.gene.list$Correlation[i],digits = 2),
                                                 " , Pvalue=",formatC(cpg.gene.list$Pvalue[i],format = "e",digits = 2),")"))+
     plot_layout(guides = "collect") & theme(legend.position = 'non')
   result.plot[[i]] <- p4
  
  }
  if(return.csv){
    return(list(cpg.gene.list=cpg.gene.list,plots=result.plot))
  }else{
    return(result.plot)
  }
}
