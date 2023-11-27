library(ggplot2)

args<-commandArgs(TRUE)
#args <- c("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/Genotyping/QC/GRCh37/QCoutput_NIH37/",
#          "NIH37")
setwd(args[1])
prefix <- args[2]

ibd <- read.table(file =paste0(prefix,"_QCd_ibd.genome"),stringsAsFactors = F,header = T )
ibd$Relatedness[ibd$RT=="UN"] <- "Unrelated"
ibd$Relatedness[ibd$RT=="HS"] <- "Half Sibs"
ibd$Relatedness[ibd$RT=="OT"] <- "Other related"
ibd$Relatedness[ibd$RT=="PO"] <- "Parent offspring"
kin <- read.table(file =paste0(prefix,"_QCd_kingX.kin0") ,stringsAsFactors = F,header = T )

outliers <- kin[kin$KinshipX >= 0.25,]

pdf(paste0(prefix,"relatedness.pdf"))

hist(kin$KinshipX,main =paste0(prefix, " KinShip Coefficient"),xlab = "Kinship Coefficient")

hist(ibd$PI_HAT,main =paste0(prefix, " IBD Pi_hat"),xlab = "Pi_hat")

ggplot(data = ibd,aes(Z0,Z1,color=Relatedness))+geom_point()+ggtitle(paste0(prefix," IBD with outliers"))
ggplot(data = ibd[!((ibd$IID1 %in% outliers$ID1)&(ibd$IID2 %in% outliers$ID2)), ],aes(Z0,Z1,color=Relatedness))+geom_point()+
  ggtitle(paste0(prefix," IBD without outliers"))
dev.off()

outliers1 <- cbind(outliers$FID1,outliers$ID1)
write.table(outliers1,file = paste0(prefix,"_RelatednessOutliers.txt"),quote = F,row.names = F,col.names = F)		
