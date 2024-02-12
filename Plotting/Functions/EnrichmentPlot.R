EnrichmentPlot <- function (result, xaxis ,yaxis, num = 5, colorby, title = "") 
{
  xaxis = match.arg(xaxis, colnames(result))
  yaxis = match.arg(yaxis, colnames(result))
  colorby = match.arg(colorby, colnames(result))
  if (!is.numeric(num)) {
    stop("num should be an integer.")
  }
  if (nrow(result) > num) {
    result = result[seq_len(num), ]
  }
  library(stringr)
  library(ggplot2)
  result$x_ <- result[,xaxis]
  result$y_ <- result[,yaxis]
  result$c_ <- as.numeric(result[,colorby])
  ggplot(result,aes(x = x_, y = y_, fill = c_)) +
    geom_bar(stat = "identity", width = 0.6) +
    ggtitle(title)+
    ylab("")+xlab(xaxis)+
    scale_colour_gradient2()+
    guides(fill=guide_legend(title = colorby))+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.1,"in"),
          title = element_text(hjust = 0,size = 8),
          aspect.ratio = 2/1)
}
