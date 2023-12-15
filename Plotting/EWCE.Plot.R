EWCE.Plot <- function(EWCE.result, q.threshold=0.05, y.lab="",x.lab="",plot.title=""){
  library(ggplot2)
  ast_q <- rep("", dim(EWCE.result)[1])
  ast_q[EWCE.result$q < q.threshold] <- "*"
  EWCE.result$ast_q <- ast_q
  EWCE.result$sd_from_mean[EWCE.result$sd_from_mean < 0] <- 0
  EWCE.result$y_ast <- EWCE.result$sd_from_mean * 1.05
  EWCE.result$abs_sd <- abs(EWCE.result$sd_from_mean)
  upperLim <- max(abs(EWCE.result$sd_from_mean), na.rm = TRUE)
  if(y.lab==""){
    y.lab="Std.Devs. from the mean"
  }
  p <- ggplot(data = EWCE.result , aes(x=CellType , y=abs_sd , fill=abs_sd)) + 
    geom_bar(stat = "identity") + 
    ylab(y.lab)+
    xlab(x.lab)+
    ggtitle(plot.title)+
    scale_fill_gradient(low = "blue", high = "red") + 
    geom_text(ggplot2::aes_string(label = "ast_q",x = "CellType", y = "y_ast"), size = 10) + 
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45,  hjust = 1)) + 
    scale_y_continuous(limits = c(-0.2 , upperLim+1))
  
  return(p)
}