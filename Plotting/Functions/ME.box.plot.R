
ME.box.plot <- function(expression1,expression2=NA, phenotype1, phenotype2=NA ,soft.power, colors, phenotype, category.col, modules, y.lab=NA, plot.title=NA, facet.col = NA,category.labels){
  
  suppressMessages(library(ggplot2))
  
  if(!identical(rownames(phenotype1) , colnames(expression1))){
    stop("Row names of phenotype1 data should be matched with colnames of expression1 data")
  }
  if(!all(is.na(expression2))){
    if(!identical(rownames(phenotype2) , colnames(expression2))){
      stop("Row names of phenotype1 data should be matched with colnames of expression1 data")
    }
    if(!all(identical(rownames(expression1) , rownames(expression2)),identical(rownames(expression1) , names(colors)),identical(rownames(expression2) , names(colors)))){
      index <- intersect(intersect(rownames(expression1),rownames(expression2)), names(colors))
      colors <- colors[index]
      expression1 <- expression1[index , ]
      expression2 <- expression2[index , ]
    }
  }

  print("Calculating ME based on expression 1...")
  MEs1 = WGCNA::moduleEigengenes(expr = t(expression1),colors = colors,softPower = soft.power)$eigengenes
  if(!all(is.na(expression2))){
  print("Calculating ME based on expression 2...")
  MEs2 = WGCNA::moduleEigengenes(expr = t(expression2),colors = colors,softPower = soft.power)$eigengenes
  }
  
  print("Generating plots...")
  results <- vector(mode = "list" , length = length(modules))
  names(results) <- modules
  for (module in modules) {
    if(!all(is.na(expression2))){
      phenotype <- data.frame(y_=c(MEs1[,paste0("ME",module)],MEs2[,paste0("ME",module)]))
      phenotype$x_ <- factor(c(phenotype1[,category.col],phenotype2[,category.col]))
    }else{
      phenotype <- data.frame(y_=MEs1[,paste0("ME",module)])
      phenotype$x_ <- factor(phenotype1[,category.col])
    }
    
    if(is.na(y.lab)){
      y.lab <- paste0(module," Module Eigengen")
    }
    
    if(!is.na(facet.col)){
      phenotype$f_ <- factor(c(phenotype1[,facet.col],phenotype2[,facet.col]), levels = c(unique(phenotype1[,facet.col]),unique(phenotype2[,facet.col])))
    }
    p <- ggplot(data = phenotype,aes(x=x_,y=y_, fill=x_))+ 
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_boxplot() + 
      scale_fill_manual(values = c("#01BEC3","#F8756D"),labels=category.labels)+
      scale_x_discrete(labels=category.labels)+
      geom_jitter(color="black", size=0.8, alpha=0.9) +
      ylab(y.lab) + 
      xlab(category.col) + 
      guides(fill=guide_legend(title=category.col)) + 
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    ylim1 = boxplot.stats(phenotype$y_)$stats[c(1, 5)]
    p = p + coord_cartesian(ylim = ylim1*1.5)
    
    if(is.na(plot.title)){
      p <- p + ggtitle(module)
    }else{
      p <- p + ggtitle(plot.title)
    }
    if(!is.na(facet.col) & !all(is.na(expression2))){
      p <- p + facet_grid(~f_)
    }
    
    results[[module]] <- p
  }
  return(results)
}
