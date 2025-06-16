##############################################################

##   Reference:  https://github.com/RitchieLab/eQTpLot
##               https://biodatamining.biomedcentral.com/articles/10.1186/s13040-021-00267-6  

##############################################################

combine_data <- function(GWAS.df, QTL.df, Genes.df, gene, trait, GWAS.SigPvalue, QTL.SigPvalue, 
                         rangebp , gbuild){
 
  QTL.data <- as.data.frame(QTL.df[which(QTL.df$gene %in% gene & QTL.df$pvalue.qtl < QTL.SigPvalue & !(is.na(QTL.df$beta.qtl)) & !(is.na(QTL.df$pvalue.qtl))), ])
  if (dim(QTL.data)[1] == 0){ 
    stop("QTL.df does not have any data for the gene ", paste(gene), " meeting your QTL.SigPvalue threshold")
  }
  startpos <- min(Genes.df[which(Genes.df$gene %in% gene), ] %>% dplyr::select(start)) - rangebp
  stoppos <- max(Genes.df[which(Genes.df$gene %in% gene ), ] %>% dplyr::select(stop)) + rangebp
  chromosome <- Genes.df$chr[Genes.df$gene %in% gene][1]
  if("phenotype" %in% names(GWAS.df)){
    GWAS.df <- GWAS.df[GWAS.df$phenotype == trait,]
  }
  
  gwas.data <- GWAS.df[which(GWAS.df$chr == chromosome  & GWAS.df$position >= startpos & GWAS.df$position <= stoppos & !(is.na(GWAS.df$pvalue.gwas)) & !(is.na(GWAS.df$beta.gwas))), ]
  
  NoFisher <- FALSE
  if (dim(gwas.data[which(gwas.data$pvalue.gwas < GWAS.SigPvalue), ])[1] == 0) {
    NoFisher <- TRUE
    print(paste("WARNING: GWAS.df does not contain any SNPs with p-value <", GWAS.SigPvalue ,"within the", 
                range, "kb flanking the gene ", gene, " for the trait ", trait, ". QTL Enrcihment Plot statistics will not be calculated"))
  }
  
  QTL.data <- dplyr::ungroup(QTL.data)
  gwas.data$snp <- as.factor(gwas.data$snp)
  QTL.data$snp <- as.factor(QTL.data$snp)
  
  Combined.Data <- dplyr::left_join(gwas.data,QTL.data, by = "snp")
  
  Combined.Data$DirectionOfEffect_GWAS <- ifelse(Combined.Data$beta.gwas < 0, "Negative", ifelse(Combined.Data$beta.gwas > 0, "Positive", NA))
  Combined.Data$DirectionOfEffect_QTL <- ifelse(Combined.Data$beta.qtl < 0, "DOWN", ifelse(Combined.Data$beta.qtl > 0, "UP", NA))
  Combined.Data$Congruence <- (Combined.Data$beta.gwas * Combined.Data$beta.qtl)
  Combined.Data$Congruence <- ifelse(Combined.Data$Congruence < 0, "Incongruent-QTL", ifelse(Combined.Data$Congruence > 0, "Congruent-QTL", NA))
  Combined.Data$NeglogQTLpValue <- -(log10(Combined.Data$pvalue.qtl))
  Combined.Data$Neglog10pvalue_GWAS <- -(log10(Combined.Data$pvalue.gwas))
  Combined.Data <- Combined.Data[which(!(is.na(Combined.Data$pvalue.gwas))), ]
  Combined.Data$significance <- ifelse(Combined.Data$pvalue.gwas >= GWAS.SigPvalue, "Non-significant", "Significant")
  Combined.Data$Congruence[is.na(Combined.Data$Congruence)] <- "Non-QTL"
  Combined.Data <- Combined.Data %>% dplyr::mutate(Congruence = factor(Congruence, levels = c("Non-QTL","Congruent-QTL", "Incongruent-QTL"), ordered = TRUE))
  Combined.Data$isQTL <- Combined.Data$Congruence
  Combined.Data$isQTL <- ifelse(Combined.Data$isQTL == "Non-QTL", "Non-QTL", "QTL")  
  
  return(list(Combined.Data=Combined.Data,Fisher=NoFisher))
}
##########################################################################################################

plot.coloc <- function(GWAS.df, QTL.df, Genes.df, gene, trait="",coloc.title="", with.legend = T ,GWAS.SigPvalue=1e-5, QTL.SigPvalue=1e-5, 
                       rangebp=5e+5 , gbuild="hg19", congruence=F, Return.CSV=F){
  suppressMessages(library(ggplot2))
  suppressMessages(library(biomaRt))
  suppressMessages(library(Gviz))
  suppressMessages(library(patchwork))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggrepel))
  suppressMessages(library(Rfast))
  
  GWAS.df$chr <- as.numeric(GWAS.df$chr)
  GWAS.df$position <- as.numeric(GWAS.df$position)
  GWAS.df$snp <- as.character(GWAS.df$snp)
  GWAS.df$pvalue <- as.numeric(GWAS.df$pvalue)
  GWAS.df$pvalue[GWAS.df$pvalue==0] <- Rfast::nth(unique(GWAS.df$pvalue),2,descending = F)/100
  GWAS.df$beta <- as.numeric(GWAS.df$beta)
  
  QTL.df$gene <- as.character(QTL.df$gene)
  QTL.df$snp <- as.character(QTL.df$snp)
  QTL.df$pvalue <- as.numeric(QTL.df$pvalue)
  QTL.df$pvalue[QTL.df$pvalue==0] <- Rfast::nth(unique(QTL.df$pvalue),2,descending = F)/100
  QTL.df$beta <- as.numeric(QTL.df$beta)
  
  Genes.df$chr <- as.numeric(Genes.df$chr)
  Genes.df$start <- as.numeric(Genes.df$start)
  Genes.df$stop <- as.numeric(Genes.df$stop)
  Genes.df$gene <- as.character(Genes.df$gene)
  
  names(GWAS.df) <- tolower(names(GWAS.df))
  names(QTL.df) <- tolower(names(QTL.df))
  names(Genes.df) <- tolower(names(Genes.df))
  
  names(GWAS.df)[names(GWAS.df)=="pvalue"] <- "pvalue.gwas"
  names(GWAS.df)[names(GWAS.df)=="beta"] <- "beta.gwas"
  
  names(QTL.df)[names(QTL.df)=="pvalue"] <- "pvalue.qtl"
  names(QTL.df)[names(QTL.df)=="beta"] <- "beta.qtl"
  
  GWAS.df <- GWAS.df[!duplicated(GWAS.df$snp),]
  
  if(!Return.CSV){
    temp <- combine_data(GWAS.df, QTL.df, Genes.df, gene, trait, GWAS.SigPvalue, QTL.SigPvalue, 
                 rangebp, gbuild)
    Combined.QTL.GWAS.Data <- temp$Combined.Data
    NoFisher <- temp$Fisher
    
    minpos <- min(Combined.QTL.GWAS.Data$position, na.rm = TRUE)
    maxpos <- max(Combined.QTL.GWAS.Data$position, na.rm = TRUE)
    miny <- min(Combined.QTL.GWAS.Data$Neglog10pvalue_GWAS , na.rm = T)
    maxy <- max(Combined.QTL.GWAS.Data$Neglog10pvalue_GWAS , na.rm = T)+1
    chromosome <- Genes.df$chr[Genes.df$gene %in% gene][1]
    
    if (dim(Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ])[1] == 0) {
      Congruentdata <- FALSE
    }else {
      Congruentdata <- TRUE
    }
    if (dim(Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ])[1] == 0) {
      Incongruentdata <- FALSE
    }else {
      Incongruentdata <- TRUE
    }
    if (gbuild == "hg19") {
      hostname <- "https://grch37.ensembl.org"
    }
    if (gbuild == "hg38") {
      hostname <- "https://apr2020.archive.ensembl.org"
    }
    bm <- biomaRt::useMart(host = hostname, biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    biomTrack <- Gviz::BiomartGeneRegionTrack(genome = gbuild, chromosome = chromosome, 
                                              start = minpos, end = maxpos, filter = list(with_refseq_mrna = TRUE), 
                                              name = "ENSEMBL", background.panel = "gray95", biomart = bm,margin = c(-3, -3))
    gtrack <- Gviz::GenomeAxisTrack(fontcolor = "#000000", fontsize = 14, margin = c(-3, -3))
    itrack <- Gviz::IdeogramTrack(genome = gbuild, chromosome = paste0("chr",chromosome),margin = c(-3, -3))
    genetracks <- patchwork::wrap_elements(panel = (grid::grid.grabExpr(Gviz::plotTracks(list(itrack,biomTrack, gtrack), collapseTranscripts = "meta", transcriptAnnotation = "symbol", 
                                                                                         chromosome = chromosome, 
                                                                                         from = min(Combined.QTL.GWAS.Data$position, na.rm = TRUE), to = max(Combined.QTL.GWAS.Data$position, na.rm = TRUE), 
                                                                                         showTitle = FALSE, labelPos = "below", distFromAxis = 10, innermargin = 0, maxHeight = (2 * 10), 
                                                                                         minHeight = (2 * 10), sizes = c(1,2, 1), margin = c(-3, -3)))))
    if(coloc.title==""){
      title1 <- paste("GWAS of ", trait, ", colored by Cis QTL data for ", paste(gene,collapse = ", "), sep = "")
    }else{
      title1 <- coloc.title
    }
    text_data <- data.frame(cpg=gene , x=Genes.df$start[Genes.df$gene %in% gene])
    text_data$y <- maxy-6
    p1 <- ggplot2::ggplot(data = Combined.QTL.GWAS.Data, aes(stroke = 0)) + 
      ggplot2::coord_cartesian(xlim = c(minpos, maxpos),ylim = c(miny,maxy), expand = FALSE) + 
      ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Non-QTL"), shape = 15, color = "black", alpha = 0.2, aes(x = position, y = Neglog10pvalue_GWAS)) + 
      ggplot2::xlab("") + ggplot2::ylab(bquote(-log10(P[GWAS]))) + 
      ggplot2::ggtitle(title1) +
      ggplot2::scale_shape_manual("GWAS Direction of Effect", values = c(Negative = 25, Positive = 24), na.value = 22) + 
      ggplot2::guides(alpha = "none", size = guide_legend("QTL Normalized Effect Size", override.aes = list(shape = 24, color = "black", fill = "grey"), title.position = "top", order = 2), 
                      shape = guide_legend(title.position = "top", direction = "vertical", order = 1, override.aes = list(size = 3, fill = "grey"))) + 
      #ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25)) + 
      ggplot2::geom_vline(xintercept = Genes.df$start[Genes.df$gene %in% gene], linetype = "twodash", color = "lightblue4", linewidth = 0.2) +
      ggplot2::geom_hline(yintercept = -log10(GWAS.SigPvalue), linetype = "solid", color = "red", linewidth = 0.5) +
      ggrepel::geom_text_repel(data = text_data,aes(x=x,y=y,label=cpg),nudge_x = -5,nudge_y = 0.9,colour="blue", angle=90) +
      ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      ggplot2::theme(plot.margin = unit(c(0, 1, -0.8, 0), "cm"),title = element_text(size = 10))
    
    if(!with.legend){
      p1 <- p1+theme(legend.position = "none")
    }
    
    if (congruence == TRUE) {
      if (Congruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Congruent-QTL"), alpha = 1, 
                                       aes(x = position, y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, 
                                           alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(beta.qtl))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[mQTL]), paste("Congruous SNPs"))), low = "#000099", 
                                       high = "#33FFFF", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.QTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogQTLpValue), na.rm = TRUE), 
                                                  max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), na.rm = TRUE)))
      }
      if (Congruentdata == TRUE & Incongruentdata == TRUE) {
        p1 <- p1 + ggnewscale::new_scale_fill()
      }
      if (Incongruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Incongruent-QTL"), 
                                       alpha = 1, aes(x = position,y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, 
                                                      alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(beta.qtl))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[mQTL]), paste("Incongruous SNPs"))), low = "#990000", 
                                       high = "#FFCC33", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.QTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogQTLpValue), na.rm = TRUE), 
                                                  max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue),na.rm = TRUE)))
      }
    }else{
      p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, isQTL == "QTL"), alpha = 1, 
                                     aes(x = position, y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(beta.qtl))) + 
        ggplot2::scale_fill_viridis_c((bquote(-log10(P[mQTL]))), option = "C", guide = guide_colorbar(title.position = "top"), 
                                      limits = c(min(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), 
                                                     na.rm = TRUE), max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), na.rm = TRUE)))
    }
    
    if (nrow(as.table(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))) < 2 | 
        ncol(as.table(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))) < 2) {
      NoFisher <- TRUE
      print("Not enough data to compute enrichment significance for Enrichment Plot")
    }
    if (NoFisher == FALSE) {
      fisher <- fisher.test(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))
      fpvalue <- fisher$p.value
    }
    
    if(congruence){
      p2 <- ggplot2::ggplot(Combined.QTL.GWAS.Data) + 
        ggplot2::aes(x = significance, y = 1, fill = Congruence) + 
        ggplot2::theme(legend.title = element_blank())+
        ggplot2::geom_bar(stat = "identity", position = "fill") + 
        ggplot2::ggtitle(paste("Enrichment of QTLs among\nGWAS-significant SNPs")) + 
        ggplot2::ylab("Proportion of SNPs\nthat are QTLs") + 
        ggplot2::xlab(paste("GWAS significance\n(threshold p <", GWAS.SigPvalue, ")")) + ggplot2::ylim(0, 1.2) + 
        if (NoFisher == FALSE) {
          ggpubr::geom_signif(y_position = c(1.1, 1.1), xmin = c("Non-significant"), xmax = c("Significant"), 
                              annotation = (paste("p =", formatC(fpvalue, format = "e", digits = 2))), tip_length = 0.05)
        }
    }else{
      p2 <- ggplot2::ggplot(Combined.QTL.GWAS.Data) + 
        ggplot2::aes(x = significance, y = 1, fill = isQTL) + 
        ggplot2::theme(legend.title = element_blank())+
        ggplot2::geom_bar(stat = "identity", position = "fill") + 
        ggplot2::ggtitle(paste("Enrichment of QTLs among\nGWAS-significant SNPs")) + 
        ggplot2::ylab("Proportion of SNPs\nthat are QTLs") + 
        ggplot2::xlab(paste("GWAS significance\n(threshold p <", GWAS.SigPvalue, ")")) + ggplot2::ylim(0, 1.2) + 
        if (NoFisher == FALSE) {
          ggpubr::geom_signif(y_position = c(1.1, 1.1), xmin = c("Non-significant"), xmax = c("Significant"), 
                              annotation = (paste("p =", formatC(fpvalue, format = "e", digits = 2))), tip_length = 0.05)
        }
    }             
    if(congruence){
      if(Congruentdata & Incongruentdata){
        
        df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
        df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
        if(nrow(df1) >= 3){
          pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation in Congruent QTLs"))
          pearson.congruent <- list(estimate = 0,p.value=1)
        }
        if(nrow(df2) >= 3){
          pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation in Incongruent QTLs"))
          pearson.Incongruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df1 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#000099") +
          ggplot2::geom_smooth( data = df1, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Congruent-QTL"), 
                                formula = y ~ x, method = "lm") +
          ggplot2::geom_point(data = df2,aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), color = "#990000") +
          ggplot2::geom_smooth(data = df2, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Incongruent-QTL"), 
                               formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#000099","#990000")) + 
          ggplot2::xlab((bquote(-log10(P[mQTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.congruent$p.value, format = "e", digits = 2)), 
                             color = "#000099", hjust = 1, vjust = 1) + 
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.Incongruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.Incongruent$p.value, format = "e", digits = 2)), 
                             color = "#990000", hjust = 1, vjust = 2.5) 
      }
      if(Congruentdata & !Incongruentdata){
        
        df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
        
        if(nrow(df1) >= 3){
          pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation in Congruent QTLs"))
          pearson.congruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df1 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#000099") +
          ggplot2::geom_smooth( data = df1, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Congruent-QTL"), 
                                formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#000099")) + 
          ggplot2::xlab((bquote(-log10(P[mQTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.congruent$p.value, format = "e", digits = 2)), 
                             color = "#000099", hjust = 1, vjust = 1)
      }
      if(!Congruentdata & Incongruentdata){
        
        df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
        
        if(nrow(df2) >= 3){
          pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation in Incongruent QTLs"))
          pearson.Incongruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df2 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#990000") +
          ggplot2::geom_point(data = df2,aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), color = "#990000") +
          ggplot2::geom_smooth(data = df2, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Incongruent-QTL"), 
                               formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#990000")) + 
          ggplot2::xlab((bquote(-log10(P[mQTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.Incongruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.Incongruent$p.value, format = "e", digits = 2)), 
                             color = "#990000", hjust = 1, vjust = 1) 
      }
    }else{
      df <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & 
                                           (Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL" | Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL")), ]
      
      if(nrow(df) >= 3){
        pearson.all <- cor.test(df$NeglogQTLpValue,df$Neglog10pvalue_GWAS, method = "pearson")
      }else{
        print(paste("Not enough data to compute pearson correlation in all QTLs"))
        pearson.all <- list(estimate = 0,p.value=1)
      }
      
      p3 <- ggplot(data =df , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color = "#360f70") +
        ggplot2::geom_smooth(data = df, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, color = Congruence), method = "lm", formula = (y ~ x), color = "#ffee00") +
        scale_colour_manual(name="Direction of effect", values=c("#ffee00")) + 
        ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.all$estimate, 3), "\np = ", formatC(pearson.all$p.value, format = "e", digits = 2)), 
                           color = "#360f70", hjust = 1, vjust = 1) 
    }
    
    p4 <- p1 + genetracks + patchwork::plot_spacer() + (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2, 3))) + 
      patchwork::plot_layout(ncol = 1, height = c(4, 2, 0.1, 2)) + 
      patchwork::plot_annotation(tag_levels = "A", tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12))
    
    p1 <- p1 + genetracks + patchwork::plot_spacer() + patchwork::plot_layout(ncol = 1, height = c(4, 2, 0.1))
    
    return(list(plot.merged = p4, plot.coloc=p1,plot.cor=p3,plot.enrich=p2)) 
    
  }else{ #Return.CSV = True
    cor.result <- as.data.frame(matrix(data = NA,nrow = length(gene), ncol = 7))
    names(cor.result) <- c("gene","congruent.cor","Congruent.Pval","Incongruent.cor","Incongruent.Pval","All.cor","All.Pval")
    for(i in 1:length(gene)){
      Combined.QTL.GWAS.Data <- combine_data(GWAS.df, QTL.df, Genes.df, gene[i], trait, GWAS.SigPvalue, QTL.SigPvalue, 
                                             rangebp, gbuild)$Combined.Data
      cor.result$gene[i] <- gene[i]
      
      df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
      df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
      df <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$beta.qtl) & 
                                          (Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL" | Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL")), ]
      if(nrow(df1) >= 3){
        pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        cor.result$congruent.cor[i] <- pearson.congruent$estimate
        cor.result$Congruent.Pval[i] <- pearson.congruent$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Congruent QTLs"))
      }
      if(nrow(df2) >= 3){
        pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        cor.result$Incongruent.cor[i] <- pearson.Incongruent$estimate
        cor.result$Incongruent.Pval[i] <- pearson.Incongruent$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Incongruent QTLs"))
      }
      if(nrow(df) >= 3){
        pearson.all <- cor.test(df$NeglogQTLpValue,df$Neglog10pvalue_GWAS, method = "pearson")
        cor.result$All.cor[i] <- pearson.all$estimate
        cor.result$All.Pval[i] <- pearson.all$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in all QTLs"))
      }
      
    }
    
    return(cor.result)
  }
  
}
