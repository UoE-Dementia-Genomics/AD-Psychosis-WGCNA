EWCE.Plot.ctd <- function (ctd, genes, level = 1, metric = "mean_exp") {
  if (tolower(metric) %in% c("expr", "exp", "expression", "mean_exp", 
                             "avgexp")) {
    metric <- "mean_exp"
    y.lab <- "Mean expression"
  }
  if (tolower(metric) == "specificity") {
    metric <- "specificity"
    y.lab <- "Specificity"
  }
  if (tolower(metric) == "specificity_quantiles") {
    metric <- "specificity_quantiles"
    y.lab <- "Quantile specificity"
  }
  metric <- stringr::str_to_sentence(metric)
  mat <- as.matrix(ctd[[level]][[tolower(metric)]])
  genes <- genes[genes %in% rownames(mat)]
  if(length(genes) > 1)
    plot_data <- reshape2::melt(mat[genes, ], id.vars = "genes")
  else{
    plot_data <- reshape2::melt(mat[genes, ], id.vars = "genes")
    plot_data$var1 <- genes
    plot_data$var2 <- row.names(plot_data)
    plot_data <- plot_data[,c(2,3,1)]
  }
  colnames(plot_data) <- c("Gene", "Celltype", metric)
  gp <- ggplot(plot_data, aes_string(x = "Celltype", y = metric, fill = metric)) + 
    scale_fill_gradient(low = "blue", high = "red",name = y.lab) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    ylab(y.lab) 
  
  if(length(genes) > 1)
    gp <- gp + facet_grid(Gene ~ .)
  else
    gp <- gp + ggtitle(genes)
  
  gp <- gp +
        theme(axis.text.x = element_text(angle = 45,  hjust = 1), strip.background = element_rect(fill = "white"), 
              strip.text = element_text(color = "black"),plot.title = element_text(hjust = 0.5))
  if (metric == "Specificity") {
    gp <- gp + scale_y_continuous(breaks = c(0, 0.5, 1), 
                                  limits = c(0, 1))
  }
  return(gp)
}