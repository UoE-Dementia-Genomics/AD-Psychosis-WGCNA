soft.power.plot <- function(soft.data, select.pow){
  
  library(ggplot2)
  
  plot.data <- data.frame(soft=soft.data$fitIndices[,1], scale= -sign(soft.data$fitIndices[,3])*soft.data$fitIndices[,2])
  p1 <- ggplot(data = plot.data , aes(x = soft , y = scale, label=soft)) + geom_point() + 
    geom_text(hjust=0.5,vjust=-0.7, color="red") + 
    ggtitle("Scale independence") + xlab("Soft Threshold (power)") +
    ylab("Scale Free Topology Model Fit,signed R^2") + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = 0.9, color="red")
  
  plot.data <- data.frame(soft = soft.data$fitIndices[,1], mean.k = soft.data$fitIndices[,5] )
  h=round(plot.data$mean.k[plot.data$soft==select.pow],digits = 0)
  p2 <- ggplot(data = plot.data , aes(x = soft , y = mean.k, label=soft)) + geom_point() + 
    geom_text(hjust=0.5,vjust=-0.7, color="red") + 
    ggtitle("Mean connectivity") + xlab("Soft Threshold (power)") +
    ylab("Mean Connectivity") + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    #geom_hline(yintercept = h , color="red", linetype="dashed")+
    scale_y_continuous(breaks=c(0,h,5000,10000, 20000, 30000))
  
  
  return(list(p.scale.free=p1 , p.mean.connect=p2))
}
