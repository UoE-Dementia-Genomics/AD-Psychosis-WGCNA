setwd("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Brain EWAS Paper")
library(vroom)
library(stringr)
library(ggplot2)
library(ggpubr)
library(png)
library(grid)
library(edgeR)
library(cowplot)
library(patchwork)
library(EWCE)

beta.pitts.regressed.file <- "./Pitts.All.Psycho.0.2.v2.lmCorrected.betas.rds"
beta.bdr.regressed.file <- "./BDR.All.Psycho.0.2.AD.braak.lmCorrected.betas.rds"
pheno.pitts.file <- "./Pitts.All.Psycho.0.2.lmCorrected.pheno.csv"
pheno.bdr.file <- "./BDR.All.Psycho.0.2.AD.wgcna.pheno.csv"
net.file <- "./Pitts.Psycho.0.2.Pow3.wgcna.network.rds"
preservation.file <- "./Pitts.BDR.All.Braak_ModulePreservation.rds"
expr.pitts.regressed.file="./Pitts.Expr.Psycho.0.2.logCPM.lmCorrected.rds"
cpg.replicated.module.file="./darkgreen.dist.1000.genes.csv"
enrichment.result.go.file="./darkgreen_GOEnrichment.csv"
enrichment.result.kegg.file="./darkgreen_KeggEnrichment.csv"
ppi.image.file="./darkgreen.1kb.genes.STRING.PPI.conf.0.7.png"
coloc.qtl.file="./darkgreen.eQTL.Cis.Pval.1e-05.Dist.1MB.ColocPlot.csv"
coloc.cpg.file="./Pitts.Psycho.0.2.CpGs.ColocPlot1.txt"
gwas.SC.file="./Results/Coloc/SC.2022.Eur.GWAS.txt"
gwas.EA.file="./EA.Okbay.2022.ColocPlot1.txt"
soft.power.data.file = "./Pitts.All.Psycho.0.2.v2.lmCorrected.wgcna.softThreshold.rds"
ewce.result.file="./darkgreen.EWCE.results.csv"
ctd.file = "./March2024/Results/EWCE/darkgreen.EWCE.ctd.rds"

gwas.pvalue  = 1e-5
qtl.pvalue = 1e-5
range = 5e+5
gbuild = "hg19"
congruence = F
soft.power <- 3

source("Scripts/Plot.Coloc.Function.R")
source("Scripts/colored.dendrogram.R")
source("Scripts/EnrichmentPlot.function.R")
source("Scripts/Methylation.Expression.CorrPlot.R")
source("Scripts/ME.box.plot.R")
source("Scripts/Merged.ModuleMembershipPlot.R")
source("Scripts/ModulePreservationPlot.R")
source("Scripts/ModuleMembershipPlot.R")
source("Scripts/WGCNA.SoftPower.Plot.R")
source("Scripts/EWCE.Plot.Enrichment.R")
source("Scripts/EWCE.Plot.ctd.R")

enrichment.go <- read.csv(enrichment.result.go.file)
names(enrichment.go) <- c("ID","Ontology","Term","Size","Overlap","Pvalue","FDR","Significant.Genes.In.Term")
enrichment.go.cc <- enrichment.go[enrichment.go$Ontology=="CC",]
enrichment.go.bp <- enrichment.go[enrichment.go$Ontology=="BP",]
enrichment.go.mf <- enrichment.go[enrichment.go$Ontology=="MF",]
enrichment.kegg <- read.csv(enrichment.result.kegg.file)
names(enrichment.kegg) <- c("ID","Description","Size","Overlap","Pvalue","FDR","Significant.Genes.In.Term")
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
cpg.replicated.module <- vroom(cpg.replicated.module.file,show_col_types = F)

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

ewce.result <- read.csv(ewce.result.file, stringsAsFactors = F, row.names = 1)
#####################################################################
#source("Scripts/CategoricalBoxPlot.R")
#merged.box.plot <- categorical.box.plot(expression = beta.merged,phenotype = pheno.merged,
#                category.col = "Psychosis",genes = names(net$colors[net$colors=="darkred"]),method = "mean",y.lab = "Average of beta values in the darkred module",
#                title = "",facet.col = "cohort",category.labels = c("AD-P","AD+P"))


merged.box.plot <- ME.box.plot(expression1 = beta.pitts.regressed,expression2 = beta.bdr.regressed,colors = net$colors,soft.power = soft.power,phenotype1 = pheno.pitts,
                               phenotype2 = pheno.bdr,category.col = "Psychosis",modules = "darkgreen",facet.col = "cohort",category.labels = c("AD-P","AD+P"))
merged.box.plot <- merged.box.plot$darkgreen
#####################################################################

merged.mm.plot <- Merged.Module.Membership.Plot(net.colors = net$colors,expr.mat1 = beta.pitts.regressed,expr.mat2 = beta.bdr.regressed,
                                                trait1 = data.frame(Psychosis=pheno.pitts$Psychosis),
                                                trait2 =data.frame(Psychosis=pheno.bdr$Psychosis),
                                                modules = "darkgreen",soft.power = 3,plot.title = "",
                                                legend.title = "Cohort",legend.value =c("PITT-ADRC","BDR") )
merged.mm.plot <- merged.mm.plot$darkgreen
#####################################################################

preserv.plot <- module.preservation.plot(MP.Data = mp, Modules = c("darkseagreen3","grey60","firebrick4","magenta","greenyellow","darkgreen"),
                                         title.medrank = "",title.zsummary = "" )
preserv.plot.med.rank <- preserv.plot$Plot.Med.Rank
preserv.plot.z.summary <- preserv.plot$Plot.Z.Summary

#####################################################################
tiff(filename = "./Fig2.WGCNA.darkgreen.tiff",units = "in",height = 10,width = 12,res = 500)
ggarrange(preserv.plot.z.summary, preserv.plot.med.rank,merged.box.plot, merged.mm.plot, nrow = 2,ncol = 2,labels = c("A","B","C","D"))
graphics.off()

####################################################################

plot.go.cc <- EnrichmentPlot(result = enrichment.go.cc,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Cellular Component")
plot.go.bp <- EnrichmentPlot(result = enrichment.go.bp,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Biological Process")
plot.go.mf <- EnrichmentPlot(result = enrichment.go.mf,num = 10,xaxis = "Overlap",yaxis = "Term",colorby = "Pvalue",title = "GO Molecular Function")
plot.kegg <- EnrichmentPlot(result = enrichment.kegg,xaxis = "Overlap",yaxis = "Description",num = 10,colorby = "Pvalue",title = "KEGG Pathway")

ppi.raster <-  rasterGrob(ppi.image, interpolate=TRUE,height = 1,width = 1)

df <- data.frame(x = c(1:10), y = c(1:10), z = paste0("A",1:10))
df$colors_ <- c("#8DD3C7","#BC80BD","#80B1D3","#FDB462","#FCCDE5","#BEBADA","#D9D9D9","#FB8072","#FFFFB3","#B3DE69")
df$labels_ <- str_wrap(enrichment.kegg$Description[1:10],width = 40)
df$CytoskapeCode <- c(1,10,5,6,8,3,9,4,2,7)

plot.ppi <- ggplot(data = df,aes(x,y,fill=z))+
  geom_point(alpha=0,pch=21,size=6)+
  guides(fill = guide_legend(override.aes = list(alpha=1)))+
  scale_fill_manual(values =df$colors_,labels=df$labels_,name="KEGG pathway")+
  annotation_custom(ppi.raster)+
  theme(panel.border = element_blank(),panel.background = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        legend.key=element_rect(fill="white",colour = NA),legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
ewc.plot <- EWCE.Plot.Enrichment(EWCE.result = ewce.result,q.threshold = 0.05)

GO.KEGG.Plot <- ggarrange(plot.go.bp,plot.go.mf,plot.go.cc,plot.kegg,labels = c("A","B","C","D"),nrow = 2,ncol = 2,common.legend = T,legend = "right")
EWC_ <- ggarrange(blankPlot,ewc.plot,blankPlot,nrow = 1,ncol = 3,labels = c("","F",""), widths = c(0.2,1,0.2))
GO.KEGG.EWC <- ggarrange(GO.KEGG.Plot,EWC_,nrow = 2, ncol = 1,heights = c(1.5,1))

tiff(filename = "./Fig3.PPI.enrichment.darkred.tif",units = "in",height = 10,width = 20,res = 500)
ggarrange(GO.KEGG.EWC,plot.ppi,labels = c("","E"),nrow = 1,ncol = 2, heights = c(1.5,2),widths = c(1.5,2))
graphics.off()

##################################################################################

result <- plot.coloc(GWAS.df = gwas.SC,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg06158994" ,trait = "Schizophrenia",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = range,gbuild = gbuild,
                     congruence = congruence,Return.CSV= F,coloc.title = "" ,with.legend = F)

plot.scz1 <- result$plot.coloc


result <- plot.coloc(GWAS.df = gwas.SC,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg26429850" ,trait = "Schizophrenia",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = range,gbuild = gbuild,
                     congruence = congruence,Return.CSV=F,coloc.title = "" ,with.legend = T)

plot.scz2 <- result$plot.coloc

result <- plot.coloc(GWAS.df = gwas.SC,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene =c("cg23071690","cg16055914") ,trait = "Schizophrenia",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = range,gbuild = gbuild,
                     congruence = congruence,Return.CSV=F,coloc.title = "" ,with.legend = F)

plot.scz3 <- result$plot.coloc

result <- plot.coloc(GWAS.df = gwas.EA,QTL.df = coloc.qtl,Genes.df = coloc.cpg,gene ="cg17252645" ,trait = "Educational Attainment",
                     GWAS.SigPvalue = gwas.pvalue,QTL.SigPvalue = qtl.pvalue,rangebp = range,gbuild = gbuild,
                     congruence = congruence,Return.CSV=F,coloc.title = "" ,with.legend = F)

plot.EA1 <- result$plot.coloc

tiff(filename = "./Fig4.darkgreen.Coloc.tif",units = "in",height = 14,width = 25,res = 500)
ggarrange(ggarrange(plot.scz1,plot.scz2,labels = c("A","B"),nrow = 1 , ncol = 2 , widths = c(0.85,1)),ggarrange(plot.scz3,plot.EA1,blankPlot,labels = c("C","D",""),nrow = 1 , ncol = 3 , widths = c(1,1,0.18)),nrow = 2,ncol = 1)
graphics.off()

###############################################################################

return.csv = F
cor.plots <- methylation.expression.corr(data.methylation = beta.pitts.regressed.matched.with.expr, data.expression = expr.pitts.regressed,
            No.Pairs.in.Plot = 3, cpg.gene.list = cpg.replicated.module, correlation.method = "pearson",
            phenotype = pheno.pitts, trait = "Psychosis", trait.labels = c("AD-P","AD+P"), return.csv = return.csv)
if(return.csv){
  results.csv <- cor.plots$cpg.gene.list
  write.csv(results.csv,file = "./darkgreen.Methyl-Expr.Corr.spearman.regressed.csv", row.names = F)
  cor.plots <- cor.plots$plots
}
ctd <- readRDS(ctd.file)
ARMC3 <- EWCE.Plot.ctd(ctd = ctd , genes = "ARMC3",metric = "mean_exp",level =1)
SPINT2 <- EWCE.Plot.ctd(ctd = ctd , genes = "SPINT2",metric = "mean_exp",level =1)
PNPLA7 <- EWCE.Plot.ctd(ctd = ctd , genes = "PNPLA7",metric = "mean_exp",level =1)

tiff(filename = "./Fig5.darkgreen.Methyl-Expr.Corr.Pearson.regressed.tif",units = "in",height = 10,width = 18,res = 500)
ggarrange(ggarrange(cor.plots[[1]]$methyl,cor.plots[[1]]$expr,cor.plots[[1]]$corr,
          cor.plots[[2]]$methyl,cor.plots[[2]]$expr,cor.plots[[2]]$corr,
          cor.plots[[3]]$methyl,cor.plots[[3]]$expr,cor.plots[[3]]$corr,
          nrow = 3,ncol = 3,common.legend = T,legend = "none",labels = c("A","D","G","B","E","H","C","F","I")),
          blankPlot,ggarrange(ARMC3, SPINT2, PNPLA7, nrow = 3,ncol = 1,common.legend = T,labels = c("J","K","L"),legend = "right"),nrow = 1,ncol = 3,widths = c(3,0.4,2))
graphics.off()

#######################################################################################
soft.power.data <- readRDS(soft.power.data.file)
colnames(beta.pitts.regressed) <- pheno.pitts$Individual_ID

p.dend <- colored.dendrogram(expression = beta.pitts.regressed)
p.soft <- soft.power.plot(soft.data = soft.power.data,select.pow = 3)
p.scale.free <- p.soft$p.scale.free
p.mean.connect <- p.soft$p.mean.connect

tiff(filename = "./SuppFig1.dendrogram.tif",units = "in",height = 15,width = 22,res = 500)
ggarrange(p.dend , ggarrange(blankPlot,p.scale.free , p.mean.connect ,blankPlot, nrow = 1 , ncol = 4, labels = c("","B","C","")) , 
          nrow = 2 , ncol = 1, heights = c(2,1),labels =c("A","") )
graphics.off()

mm.plots <- Module.Membership.Plot(net.colors = net$colors,expr.mat = beta.pitts.regressed,trait = data.frame(Psychosis=pheno.pitts$Psychosis),
                       modules = c("darkseagreen3","grey60","firebrick4","magenta","greenyellow","darkgreen"),size.threshold =10000,soft.power = 3,plot.title = "")

merged.box.plot <- ME.box.plot(expression1 = beta.pitts.regressed,colors = net$colors,soft.power = soft.power,phenotype1 = pheno.pitts,plot.title = "",y.lab = "Module Eigengene",
                               category.col = "Psychosis",modules = c("darkseagreen3","grey60","firebrick4","magenta","greenyellow","darkgreen"),category.labels = c("AD-P","AD+P"))

p1 <- merged.box.plot$darkseagreen3 + mm.plots$Plots$darkseagreen3+ plot_annotation(title = "darkseagreen3")+
      plot_layout(guides = "collect")& theme(legend.position = 'non')
p2 <- merged.box.plot$grey60 + mm.plots$Plots$grey60+ plot_annotation(title = "grey60")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p3 <- merged.box.plot$firebrick4 + mm.plots$Plots$firebrick4+ plot_annotation(title = "firebrick4")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p4 <- merged.box.plot$magenta + mm.plots$Plots$magenta+ plot_annotation(title = "magenta")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p5 <- merged.box.plot$greenyellow + mm.plots$Plots$greenyellow+ plot_annotation(title = "greenyellow")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')
p6 <- merged.box.plot$darkgreen + mm.plots$Plots$darkgreen+ plot_annotation(title = "darkgreen")+
  plot_layout(guides = "collect")& theme(legend.position = 'non')

tiff(filename = "./SuppFig2.other.modules.tif",units = "in",height = 14,width = 18,res = 500)
ggarrange(p1,p2,p3,p4,p5,p6,nrow = 3,ncol = 2,common.legend = T)
graphics.off()
