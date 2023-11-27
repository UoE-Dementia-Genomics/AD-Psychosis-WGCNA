## script to plot PCA output and colour by 1000G ethnicity
library(ggplot2)


calcPopDist<-function(dat.pca, ref){
# dat.pca is test individual's data across pcs
# ref is matrix with 1 row per population and 1 column per pca
	popDist<-rep(NA, nrow(ref))
	names(popDist)<-rownames(ref)
	for(i in 1:nrow(ref)){
		diffs<-dat.pca - ref[i,]
		sqdiffs<-diffs^2
		popDist[i]<-sqrt(sum(sqdiffs))
	}
	return(popDist)
}


args<-commandArgs(TRUE)

#args <- c("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/Genotyping/QC/GRCh37",
#          "NIH37","C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/Genotyping/Ref/GRCh37_HG19_Referance")

setwd(args[1])
prefix<-args[2]
refFolder<-args[3]

pcas<-read.table(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_mergedw1000G.pca.eigenvec")), stringsAsFactors = FALSE)
pop <- read.csv(paste0(refFolder,'/phase3_sample_info.csv'), header = T, stringsAsFactors = F)

pcs <- pcas
pcs <- pcs[,c(-1)]
colnames(pcs)[2:dim(pcs)[2]] <- paste0("PC",c(1:(dim(pcs)[2]-1)))
colnames(pcs)[1] <- "sample"
pcs$population <- pop$Super_Population[match(pcs$sample, pop$Sample)]
pcs[is.na(pcs)] <- " "
pcs$population[which(pcs$population== " ")]<-"NIH_Samples"
pcs$population <- as.factor(pcs$population)
pcs_our <- pcs[pcs$population=="NIH_Samples",]

pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_PCAplot_All.pdf")), width = 10, height = 10)
ggplot(pcs, aes(PC1, PC2, color= population),label=sample)+
  geom_point() + 
  scale_color_manual(values = c("AFR" = "brown1", "AMR" = "darkgoldenrod1","EAS"="darkolivegreen1","SAS"="darkorchid1","NIH_Samples"="blue","EUR"="cyan")) +
  geom_point(data=pcs_our[pcs_our$PC2< 0.015, ],pch=21, fill=NA, size=4, colour="blue", stroke=1)

ggplot(pcs, aes(PC2, PC3, color= population),label=sample)+
  geom_point() + 
  scale_color_manual(values = c("AFR" = "brown1", "AMR" = "darkgoldenrod1","EAS"="darkolivegreen1","SAS"="darkorchid1","NIH_Samples"="blue","EUR"="cyan")) +
  geom_point(data=pcs_our[pcs_our$PC2< 0.015, ],pch=21, fill=NA, size=4, colour="blue", stroke=1)

dev.off()



pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_PCAplot.pdf")), width = 10, height = 10)
ggplot(pcs_our, aes(PC1, PC2, color= population),label=sample)+
  geom_point() +
  geom_point(data=pcs_our[pcs_our$PC2< 0.015, ],pch=21, fill=NA, size=4, colour="blue", stroke=1)
ggplot(pcs_our, aes(PC2, PC3, color= population),label=sample)+
  geom_point() +
  geom_point(data=pcs_our[pcs_our$PC2< 0.015, ],pch=21, fill=NA, size=4, colour="blue", stroke=1)
dev.off()

KGped<-read.table(paste(refFolder, "/20130606_g1k.ped", sep = ""), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
popInfo<-read.table(paste(refFolder, "/PopInfo.txt", sep = ""), stringsAsFactors = FALSE, header = TRUE, sep = "\t") ## table made from 1000G website

KGped<-KGped[match(pcas[,2], KGped[,2]),]
nPops<-length(table(KGped$Population))
popInfo<-popInfo[match(popInfo$Population.Code,levels(as.factor(KGped$Population))),]

KGped<-cbind(KGped,popInfo$Super.Population.Code[match(KGped$Population, popInfo$Population.Code)])
colnames(KGped)[ncol(KGped)]<-"SuperPopulation"
nSuperPops<-length(table(KGped$SuperPopulation))



ptType<-c(20,3)[as.factor(is.na(KGped[,1]))]
ptCol<-rainbow(nPops)[as.factor(KGped$Population)]
ptCol[is.na(ptCol)]<-"black"

pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_PCAplotwith1KG.pdf")), width = 10, height = 10)
par(mfrow = c(2,2))
par(mar = c(4,4,0.75,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nPops), levels(as.factor(KGped$Population)), cex = 1.5, ncol = 3)


## alternatively plot "super populations"


ptCol<-rainbow(nSuperPops)[as.factor(KGped$SuperPopulation)]
ptCol[is.na(ptCol)]<-"black"
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)

dev.off()



## for each super population calculate cluster medians
print("Calculating population means")
nMatches<-rep(NA,20) 
for(nPCs in 2:20){
	pop.medians<-apply(pcas[,-c(1:2)][,1:nPCs], 2,aggregate, by = list(KGped$SuperPopulation), median)
	pop.medians<-cbind.data.frame(pop.medians)
	rownames(pop.medians)<-pop.medians[,1]
	pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

	## for each individual compare to each super population and find most similar
	popDistsAll<-apply(pcas[,-c(1:2)][,1:nPCs], 1, calcPopDist, pop.medians)
	popDistsAll<-t(popDistsAll)
	predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]

	compTrue<-table(predPop, KGped$SuperPopulation)
	nMatches[nPCs]<-sum(diag(compTrue))
  print(nPCs)
}


pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_SelectOptimalnPCsForPopulationPrediction.pdf")))
plot(1:20,nMatches/sum(!is.na(KGped$SuperPopulation))*100, xlab = "nPCs", ylab = "Percentage Correct")
dev.off()

nPCs<-which.max(nMatches)
pop.medians<-apply(pcas[,-c(1:2)][,1:nPCs], 2,aggregate, by = list(KGped$SuperPopulation), median)
pop.medians<-cbind.data.frame(pop.medians)
rownames(pop.medians)<-pop.medians[,1]
pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

## for each individual compare to each super population and find most similar
popDistsAll<-apply(pcas[,-c(1:2)][,1:nPCs], 1, calcPopDist, pop.medians)
popDistsAll<-t(popDistsAll)
predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]

## calculate a quality score for prediction
## ideally want one population much closer than the others
rangeDist<-t(diff(apply(popDistsAll,1,range)))
qsPred<-(apply(popDistsAll,1,quantile, 0.25)-apply(popDistsAll,1,min))/rangeDist
pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_BoxplotPrePopQCScores.pdf")), width = 12, height = 6)
par(mfrow = c(1,2))
boxplot(qsPred ~ KGped$SuperPopulation, col = rainbow(5), xlab = "Known populations")
boxplot(qsPred ~ predPop, col = rainbow(5), xlab = "Predicted populations")
dev.off()

## can define thresholds for each population based on 1000 genomes samples
#pop99<-aggregate(popDistsAll, by = list(KGped$SuperPopulation), quantile, 0.95)
#pop99Thres<-diag(as.matrix(pop99[,-1]))
#names(pop99Thres)<-colnames(popDistsAll)
#table(apply(popDistsAll, 1, min) < pop99Thres[predPop],as.factor(predPop))

print("Plotting predicted populations")

ptCol<-rainbow(nSuperPops)[as.factor(predPop)]
pdf(paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_PCAplotwith1KGpredictedPopulations.pdf")), width = 10, height = 10)
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 20),3], pcas[which(ptType == 20),4], pch = 20, col = ptCol[which(ptType == 20)])
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 3),3], pcas[which(ptType == 3),4], pch = 3, col = ptCol[which(ptType == 3)])

dev.off()


## look at predications of our sample
outPred<-cbind(pcas[,1:2], predPop, qsPred, popDistsAll)
outPred<-outPred[which(ptType == 3),]
write.csv(table(outPred$predPop), paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_TablePredictedPopulations.csv")))
write.csv(outPred, paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_PredictedPopulations.csv")), quote = FALSE, row.names = FALSE)

fam <- read.table(paste0(paste0("QCoutput_",prefix),paste0("/",prefix,"_QCd.fam")), stringsAsFactors = FALSE)

outliers <- pcs_our$sample[pcs_our$PC2< 0.015]
outliers1 <- fam[fam$V2 %in% outliers,]
write.table(outliers1[,c(1,2)],file = paste0(paste0("merge1KG_",prefix),paste0("/",prefix,"_EthnicityOutliers.txt")),quote = F,row.names = F,col.names = F)		
