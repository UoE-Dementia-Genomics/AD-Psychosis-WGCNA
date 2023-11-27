# Ref: https://r-graph-gallery.com/31-custom-colors-in-dendrogram.html

colored.dendrogram <- function(expression, phenotype=NA, color.column=NA, plot.title="", legend.title="", Colors=NA, 
                               legend.x=30, legend.y=30, plot.font.size=0.8, legend.font.size=1){
  if(!is.na(color.column)){
    if(!identical(colnames(expression) , rownames(phenotype))){
      stop("Colnames of expression should be matched with rownames of phenotype")
    }
  }
  print("Calculating distance between samples...")
  distance <- dist(t(expression), method = "euclidean")
  print("Generating hclust...")
  hc <- hclust(distance,method = "average")
  dhc <- as.dendrogram(hc,hang = 0.2)
  print("Generating plot...")
  if(!is.na(color.column)){
    colors_ <- data.frame(label=unique(phenotype[,color.column]),color=rainbow(length(unique(phenotype[,color.column]))))
    if(!is.na(Colors)[1]){
      colors_$color = Colors
    }
  }

  if(legend.title==""){
    if(is.numeric(color.column)){
      legend.title <- colnames(phenotype)[color.column]
    }else{
      legend.title <- color.column
    }
  }
  colLab<<-function(n){
    if(is.leaf(n)){
      
      #I take the current attributes
      a=attributes(n)
      
      #I deduce the line in the original data, and so the treatment and the specie.
      ligne=match(attributes(n)$label,colnames(expression))
      treatment=phenotype[ligne,color.column];
      
      for(i in 1:nrow(colors_))
      if(treatment==colors_$label[i]){
        col_treatment=colors_$color[i]
        }
      
      #Modification of leaf attribute
      attr(n,"nodePar")<-c(a$nodePar,list(pch=20,col=col_treatment,lab.col=col_treatment))
    }
    return(n)
  }
  
  # Finally I just have to apply this to my dendrogram
  if(!is.na(color.column)){
    dL <- dendrapply(dhc, colLab)
    # And the plot
    par(cex=plot.font.size)
    plot(dL , main=plot.title, ylab="Height",xlab = "", sub = "")
    par(cex=legend.font.size)
    legend(legend.x,legend.y, title = legend.title,
           legend = colors_$label, 
           col =colors_$color , 
           pt.cex = 1.5, bty = "n", pch = rep(20 , times=nrow(colors_)),
           text.col = "black", horiz = FALSE, inset = c(0, 0.1))
  }else{
    par(cex=plot.font.size)
    plot(hc, main=plot.title , xlab = "", sub = "")
  }

  
  p <- recordPlot()
  
  return(p)
  
}

