library(vroom)
library(stringr)
library(ggplot2)
library(ggpubr)
library(png)
library(grid)
library(edgeR)
library(cowplot)
library(patchwork)

args <- commandArgs(T)

beta.pitts.regressed.file <- args[1]
beta.bdr.regressed.file <- args[2]
pheno.pitts.file <- args[3]
pheno.bdr.file <- args[4]
net.file <- args[5]
preservation.file <- args[6]
expr.pitts.regressed.file <- args[7]
cpg.darkred.file <- args[8]
enrichment.result.go.file <- args[9]
enrichment.result.kegg.file <- args[10]
ppi.image.file <- args[11]
coloc.qtl.file <- args[12]
coloc.cpg.file <- args[13]
gwas.SC.file <- args[14]
gwas.EA.file <- args[16]
soft.power.data.file <- args[17]

gwas.pvalue  <- args[18]
qtl.pvalue <- args[19]
soft.power <- args[20]

source("Plotting/Plot.Coloc.Function.R")
source("Plotting/colored.dendrogram.R")
source("Plotting/EnrichmentPlot.function.R")
source("Plotting/Methylation.Expression.CorrPlot.R")
source("Plotting/ME.box.plot.R")
source("Plotting/Merged.ModuleMembershipPlot.R")
source("Plotting/ModulePreservationPlot.R")
source("Plotting/ModuleMembershipPlot.R")
source("Plotting/WGCNA.SoftPower.Plot.R")

enrichment.go <- read.csv(enrichment.result.go.file) # names(enrichment.go) <- c("ID","Ontology","Term","Size","Overlap","Pvalue","FDR","Significant.Genes.In.Term")
enrichment.go.cc <- enrichment.go[enrichment.go$Ontology=="CC",]
enrichment.go.bp <- enrichment.go[enrichment.go$Ontology=="BP",]
enrichment.go.mf <- enrichment.go[enrichment.go$Ontology=="MF",]
enrichment.kegg <- read.csv(enrichment.result.kegg.file) #names(enrichment.kegg) <- c("ID","Description","Size","Overlap","Pvalue","FDR","Significant.Genes.In.Term")
ppi.image <- readPNG(source = ppi.image.file)

coloc.qtl <- vroom(file = coloc.qtl.file,col_types = list(snp='c'))
coloc.cpg <- vroom(coloc.cpg.file,show_col_types = F)
coloc.cpg = coloc.cpg[coloc.cpg$gene %in% coloc.qtl$gene,]
gwas.SC <- vroom(gwas.SC.file,col_types = list(snp='c'))
gwas.EA <- vroom(gwas.EA.file,col_types = list(snp='c'))

beta.pitts.regressed <- readRDS(beta.pitts.regressed.file)
beta.bdr.regressed <- readRDS(beta.bdr.regressed.file)
pheno.pitts <- read.csv(pheno.pitts.file, row.names = 1 , stringsAsFactors = F)
pheno.bdr <- read.csv(pheno.bdr.file, row.names = 1 , stringsAsFactors = F)
net <- readRDS(net.file)
mp <- readRDS(preservation.file)
expr.pitts.regressed <- readRDS(expr.pitts.regressed.file)
cpg.darkred <- vroom(cpg.darkred.file,show_col_types = F)

identical(colnames(beta.pitts.regressed) , rownames(pheno.pitts))

index <- match(colnames(expr.pitts.regressed),pheno.pitts$Individual_ID)
colnames(expr.pitts.regressed) <- pheno.pitts$Basename[index]
index <- intersect(colnames(expr.pitts.regressed),colnames(beta.pitts.regressed))
expr.pitts.regressed <- expr.pitts.regressed[,index]
beta.pitts.regressed.matched.with.expr <- beta.pitts.regressed[,index]

identical(colnames(beta.pitts.regressed.matched.with.expr) , colnames(expr.pitts.regressed))
identical(colnames(beta.bdr.regressed) , rownames(pheno.bdr))
identical(colnames(beta.pitts.regressed) , rownames(pheno.pitts))

pheno.bdr$cohort = "BDR"
pheno.pitts$cohort = "PITT-ADRC"
index <- intersect(names(pheno.bdr) , names(pheno.pitts))
pheno.merged <- rbind.data.frame(pheno.bdr[,index],pheno.pitts[,index])

index <- intersect(rownames(beta.bdr.regressed) , rownames(beta.pitts.regressed))
beta.merged <- cbind.data.frame(beta.bdr.regressed[index, ],beta.pitts.regressed[index, ])

identical(colnames(beta.merged) , rownames(pheno.merged))

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

#####################################################################
merged.box.plot <- ME.box.plot(expression1 = beta.pitts.regressed,expression2 = beta.bdr.regressed,colors = net$colors,soft.power = soft.power,phenotype1 = pheno.pitts,
                               phenotype2 = pheno.bdr,category.col = "Psychosis",modules = "darkred",facet.col = "cohort",category.labels = c("AD-P","AD+P"))
#####################################################################

merged.mm.plot <- Merged.Module.Membership.Plot(net.colors = net$colors,expr.mat1 = beta.pitts.regressed,expr.mat2 = beta.bdr.regressed,
                                                trait1 = data.frame(Psychosis=pheno.pitts$Psychosis),
                                                trait2 =data.frame(Psychosis=pheno.bdr$Psychosis),
                                                modules = "darkred",soft.power = 3,plot.title = "",
                                                legend.title = "Cohort",legend.value =c("PITT-ADRC","BDR") )
merged.mm.plot <- merged.mm.plot$darkred
#####################################################################

preserv.plot <- module.preservation.plot(MP.Data = mp, Modules = c("mediumpurple3","firebrick3",
                                                                   "darkgoldenrod1",
                                                                   "darkred",
                                                                   "lightblue1"),title.medrank = "",title.zsummary = "" )
preserv.plot.med.rank <- preserv.plot$Plot.Med.Rank
preserv.plot.z.summary <- preserv.plot$Plot.Z.Summary

#####################################################################
tiff(filename = "Fig2_WGCNA.darkred1.tif",units = "in",height = 10,width = 12,res = 500)
ggarrange(merged.box.plot, merged.mm.plot, preserv.plot.z.summary, 
          preserv.plot.med.rank,nrow = 2,ncol = 2,labels = c("A","B","C","D"))

graphics.off()

####################################################################

plot.go.cc <- EnrichmentPlot(result = enrichment.go.cc,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Cellular Component")
plot.go.bp <- EnrichmentPlot(result = enrichment.go.bp,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Biological Process")
plot.go.mf <- EnrichmentPlot(result = enrichment.go.mf,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Molecular Function")
plot.kegg <- EnrichmentPlot(result = enrichment.kegg,xaxis = "Overlap",yaxis = "Description",num = 10,colorby = "Pvalue",title = "KEGG Pathway")

ppi.raster <-  rasterGrob(ppi.image, interpolate=TRUE,height = 1,width = 1)

df <- data.frame(x = c(1:10), y = c(1:10), z = paste0("A",1:10))
colors_ <- c("#8DD3C7","#BC80BD","#80B1D3","#FCCDE5","#BEBADA","#D9D9D9","#FDB462","#FB8072","#FFFFB3","#B3DE69")
labels_ <- str_wrap(enrichment.kegg$Description[1:10],width = 40)
plot.ppi <- ggplot(data = df,aes(x,y,fill=z))+
  geom_point(alpha=0,pch=21,size=6)+
  guides(fill = guide_legend(override.aes = list(alpha=1)))+
  scale_fill_manual(values =colors_,labels=labels_,name="KEGG pathway")+
  annotation_custom(ppi.raster)+
  theme(panel.border = element_blank(),panel.background = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        legend.key=element_rect(fill="white"),legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))

tiff(filename = "Fig3_PPI.enrichment.darkred.tif",units = "in",height = 8,width = 18,res = 500)
ggarrange(ggarrange(plot.kegg,plot.go.bp,plot.go.mf,plot.go.cc,labels = c("A","B","C","D"),nrow = 2,ncol = 2,common.legend = T,legend = "right"),
          plot.ppi,labels = c("","E"),nrow = 1,ncol = 2, heights = c(1.5,2),widths = c(1.5,2))
graphics.off()

##################################################################################

result <- plot.coloc(GWAS.df = gwas.SC,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg06158994" ,trait = "Schizophrenia",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = 5e+5,gbuild = "hg19",
                     congruence = congruence,Return.CSV= F)

plot.scz1 <- result$plot.coloc


result <- plot.coloc(GWAS.df = gwas.SC,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg19882179" ,trait = "Schizophrenia",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = 5e+5,gbuild = "hg19",
                     congruence = congruence,Return.CSV=F)

plot.scz2 <- result$plot.coloc


result <- plot.coloc(GWAS.df = gwas.EA,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg17252645" ,trait = "Educational Attainment",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = 5e+5,gbuild = "hg19",
                     congruence = congruence,Return.CSV=F)

plot.EA1 <- result$plot.coloc

tiff(filename = "Fig4.darkred.Coloc.tif",units = "in",height = 14,width = 20,res = 500)
ggarrange(ggarrange(plot.scz1,plot.scz2,nrow = 1,ncol = 2,labels = c("A","B")),
          ggarrange(blankPlot,plot.EA1,blankPlot,nrow = 1,ncol = 3,widths = c(1,3.2,0.75),labels = c("","C","")),
          nrow = 2,ncol = 1)
graphics.off()

###############################################################################

return.csv = F
cor.plots <- methylation.expression.corr(data.methylation = beta.pitts.regressed.matched.with.expr, data.expression = expr.pitts.regressed,
            No.Pairs.in.Plot = 2, cpg.gene.list = cpg.darkred, correlation.method = "pearson",
            phenotype = pheno.pitts, trait = "Psychosis", trait.labels = c("AD-P","AD+P"), return.csv = return.csv)
if(return.csv){
  results.csv <- cor.plots$cpg.gene.list
  write.csv(results.csv,file = "Fig5_darkred.Methyl-Expr.Corr.spearman.regressed.csv", row.names = F)
  cor.plots <- cor.plots$plots
}

tiff(filename = "Fig5_darkred.Methyl-Expr.Corr.Pearson.regressed.tif",units = "in",height = 8,width = 10,res = 500)
ggarrange(cor.plots[[1]],cor.plots[[2]],nrow = 2,ncol = 1)
graphics.off()

#######################################################################################
soft.power.data <- readRDS(soft.power.data.file)
colnames(beta.pitts.regressed) <- pheno.pitts$Individual_ID
p.dend <- colored.dendrogram(expression = beta.pitts.regressed)
p.soft <- soft.power.plot(soft.data = soft.power.data,select.pow = 3)
p.scale.free <- p.soft$p.scale.free
p.mean.connect <- p.soft$p.mean.connect

tiff(filename = "SuppFig1.dendrogram.tif",units = "in",height = 15,width = 22,res = 500)
ggarrange(p.dend , ggarrange(blankPlot,p.scale.free , p.mean.connect ,blankPlot, nrow = 1 , ncol = 4, labels = c("","B","C","")) , 
          nrow = 2 , ncol = 1, heights = c(2,1),labels =c("A","") )
graphics.off()

mm.plots <- Module.Membership.Plot(net.colors = net$colors,expr.mat = beta.pitts.regressed,trait = data.frame(Psychosis=pheno.pitts$Psychosis),
                       modules = c("firebrick3","mediumpurple3","lightblue1","darkgoldenrod1","darkred"),size.threshold = 1500,soft.power = 3,return.members = F,plot.title = "")


merged.box.plot <- ME.box.plot(expression1 = beta.pitts.regressed,colors = net$colors,soft.power = soft.power,phenotype1 = pheno.pitts,plot.title = "",y.lab = "Module Eigengene",
                               category.col = "Psychosis",modules = c("firebrick3","mediumpurple3","lightblue1","darkgoldenrod1","darkred"),category.labels = c("AD-P","AD+P"))

p1 <- merged.box.plot$firebrick3 + mm.plots$firebrick3+ plot_annotation(title = "firebrick3")+
      plot_layout(guides = "collect")& theme(legend.position = 'non')
p2 <- merged.box.plot$mediumpurple3 + mm.plots$mediumpurple3+ plot_annotation(title = "mediumpurple3")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p3 <- merged.box.plot$lightblue1 + mm.plots$lightblue1+ plot_annotation(title = "lightblue1")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p4 <- merged.box.plot$darkgoldenrod1 + mm.plots$darkgoldenrod1+ plot_annotation(title = "darkgoldenrod1")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p5 <- merged.box.plot$darkred + mm.plots$darkred+ plot_annotation(title = "darkred")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')

tiff(filename = "SuppFig2.other.modules.tif",units = "in",height = 14,width = 18,res = 500)
ggarrange(p1,p2,p3,p4,p5,nrow = 3,ncol = 2,common.legend = T)
graphics.off()


