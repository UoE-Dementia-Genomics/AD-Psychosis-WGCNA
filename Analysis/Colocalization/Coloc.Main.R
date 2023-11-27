
get.lead.snps <- function(loci.data){
  result <- vector(mode="character" , length=0)
  for (i in 1:nrow(loci.data)) {
    if((as.numeric(loci.data[i,]$END) - as.numeric(loci.data[i,]$START)) < 2){
      result <- append(result,loci.data[i,]$SNP)# Lead SNP
    }
  }
  return(data.frame(id=result))
}
##################################################################
args<-commandArgs(TRUE)

QTL.file <- args[1]
GWAS.file <- args[2]
loci.file <- args[3]
ref.genome.prefix <- args[4]
distance_ <- as.numeric(args[5])
type_ <- args[6]
ld.threshold <- as.numeric(args[7])
out.pref <- args[8]
#################################################################
print("Input arguments:")
print(paste0("     QTL file= ",QTL.file))
print(paste0("     GWAS Summary statistic file= ",GWAS.file))
print(paste0("     Loci file= ",loci.file))
print(paste0("     Coloc type= ",type_))
print(paste0("     Distance for selecting region based on LD (if the region length < 2)= ",distance_))
print(paste0("     Reference genome binary files name (--bfile option in plink)= ",ref.genome.prefix))
print(paste0("     Output prefix= ",out.pref))
cat('\n')

print("Loading libraries...")
suppressMessages(library(coloc))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
options(datatable.fread.datatable=FALSE)

print("Reading QTL file...")
QTLs <- fread(QTL.file, stringsAsFactors = F, header = T)
names(QTLs) <- toupper(names(QTLs))
QTLs <- QTLs[!duplicated(QTLs$SNP),] 
QTLs$MAF <- as.numeric(QTLs$MAF)
QTLs$SE <- as.numeric(QTLs$SE)
QTLs$BETA <- as.numeric(QTLs$BETA)
QTLs$N <- as.numeric(QTLs$N)
QTLs <- QTLs[(!is.na(QTLs$BETA))&(!is.na(QTLs$SE))&(!is.na(QTLs$SNP))&(!is.na(QTLs$MAF)),]
QTLs <- QTLs[QTLs$MAF>0,]
QTLs <- QTLs[QTLs$MAF<1,]
QTLs_coloc = list(beta = as.numeric(QTLs$BETA), # Beta values for allele1 of each CPG
                  varbeta = QTLs$SE^2, # Variance of each beta value
                  type = type_, # quant for quantitative or cc for binary outcome
                  snp = QTLs$SNP, # IDs of QTLs. MUST be the same as in the GWAS dataset
                  MAF = QTLs$MAF, # Frequency of allele1 of each CPG
                  N = QTLs$N) # Sample size of the associated study.In general will be the max of available sample size values
check_ <- check_dataset(QTLs_coloc) # That functions returns NULL if everything is OK for the coloc analysis
if(is.null(check_)){
  print("QTL dataset is OK")
}

print("Reading GWAS summary statistics...")
GWAS <- fread(file=GWAS.file,stringsAsFactors=F,header = T)
names(GWAS) <- toupper(names(GWAS))
GWAS <- GWAS[!duplicated(GWAS$SNP),]
GWAS$MAF <- as.numeric(GWAS$MAF)
GWAS$SE <- as.numeric(GWAS$SE)
GWAS$BETA <- as.numeric(GWAS$BETA)
GWAS$N <- as.numeric(GWAS$N)
GWAS <- GWAS[(!is.na(GWAS$BETA))&(!is.na(GWAS$SE))&(!is.na(GWAS$SNP))&(!is.na(GWAS$MAF)),]
GWAS <- GWAS[GWAS$MAF>0,]
GWAS <- GWAS[GWAS$MAF<1,]
GWAS_coloc = list(beta = GWAS$BETA, 
                  varbeta = GWAS$SE^2,
                  type = type_,
                  snp = GWAS$SNP, 
                  MAF = GWAS$MAF,
                  N = GWAS$N)
check_ <- check_dataset(GWAS_coloc)
if(is.null(check_)){
  print("GWAS dataset is OK")
}

print("Reading loci file...")
loci <- fread(file = loci.file,stringsAsFactors = F,header = T)
names(loci) <- toupper(names(loci))
loci$START <- as.numeric(loci$START)
loci$END <- as.numeric(loci$END)

snp.list <- get.lead.snps(loci)
write.table(snp.list,file = paste0(loci.file,"lead.snp.list"),row.names = F,col.names = F,quote = F)
ld <- NULL
if(nrow(snp.list) > 0){
  print("Calculating LD...")
  plink_cmd <- paste("plink --r2 --bfile",ref.genome.prefix , "--ld-snp-list", paste0(loci.file,"lead.snp.list") ,"--ld-window-r2",ld.threshold,"--ld-window-kb",distance_,"--out",loci.file)
  if(file.exists(paste0(loci.file,".ld"))){
    ld <- fread(paste0(loci.file,".ld"),stringsAsFactors = F , header = T)
  }else{
    exit_code <- system(plink_cmd, ignore.stdout=T,wait = T)
    rm <- file.remove(paste0(loci.file,".log"))
    if(!(exit_code>0)){
      ld <- fread(paste0(loci.file,".ld"),stringsAsFactors = F , header = T)
    }
  } 
}
rm <- file.remove(paste0(loci.file,"lead.snp.list"))

log_ <- as.data.frame(matrix(data = "",nrow = nrow(loci),ncol = 23))
names(log_) <- c("Coloc.Type","LD.Threshold","Distance","UniqCpG.QTL","UniqSNPs.QTL","UniqSNPs.GWAS","CommonSNPs.QTL.GWAS",
                 "locus","locus.Start","locus.End","locus.length","UniqSNPs.locus","CommonSNPs.QTL.locus","PP.H0.abf","PP.H1.abf",
                 "PP.H2.abf","PP.H3.abf","PP.H4.abf","sum.H3.H4","SNPs","nSNPs","CpGs","nCpGs")

result <- vector(mode = "list",length = nrow(loci))
names(result) <- paste0("locus.",loci$SNP)
for (i in 1:nrow(loci)) {
  print(paste("Processing locus",loci$SNP[i]))
  log_$nSNPs[i] = 0
  log_$PP.H0.abf[i] = 0
  log_$PP.H1.abf[i] = 0
  log_$PP.H2.abf[i] = 0
  log_$PP.H3.abf[i] = 0
  log_$PP.H4.abf[i] = 0
  log_$sum.H3.H4[i] = 0
  log_$locus[i] = loci$SNP[i]
  if(!is.null(ld)){
    if(loci$SNP[i] %in% ld$SNP_A){
      ld1 <- ld[ld$SNP_A==loci$SNP[i],]
      loci$START[i] <- min(ld1$BP_B)
      loci$END[i] <- max(ld1$BP_B)
    }
  }
  log_$locus.Start[i] <- loci$START[i]
  log_$locus.End[i] <- loci$END[i]
  log_$locus.length[i] <- loci$END[i] - loci$START[i]
  index <- (GWAS$CHR == loci$CHR[i]) & (as.numeric(GWAS$POS) > loci$START[i]) & (as.numeric(GWAS$POS) < loci$END[i])
  GWAS.loci <- GWAS[index , ]
  c <- length(QTLs$SNP[QTLs$SNP %in% GWAS.loci$SNP])
  log_$CommonSNPs.QTL.locus[i] <- c
  log_$UniqSNPs.locus[i] <- length(unique(GWAS.loci$SNP))
  if(c > 0){
    GWAS_coloc = list(beta = GWAS.loci$BETA,
                      varbeta = GWAS.loci$SE^2,
                      type = type_,
                      snp = GWAS.loci$SNP, 
                      MAF = GWAS.loci$MAF,
                      N = GWAS.loci$N)
    
    result.coloc = coloc.abf(GWAS_coloc, QTLs_coloc)
    result.QTL = QTLs[QTLs$SNP %in% GWAS.loci$SNP,]
    result[[i]] <- list(result.coloc = result.coloc,result.QTL = result.QTL)
    temp = as.data.frame(result.coloc$summary)
    log_$nSNPs[i] = temp["nsnps",1]
    log_$PP.H0.abf[i] = temp["PP.H0.abf",1]
    log_$PP.H1.abf[i] = temp["PP.H1.abf",1]
    log_$PP.H2.abf[i] = temp["PP.H2.abf",1]
    log_$PP.H3.abf[i] = temp["PP.H3.abf",1]
    log_$PP.H4.abf[i] = temp["PP.H4.abf",1]
    log_$sum.H3.H4[i] = temp["PP.H3.abf",1] + temp["PP.H4.abf",1]
    cpgs = paste(unique(result.QTL$CPG),collapse=',')
    log_$CpGs[i] = cpgs
    log_$nCpGs[i] = length(unique(result.QTL$CPG))
    log_$SNPs[i] <- paste(unique(result.coloc$results$snp),collapse=',')
    log_$nSNPs[i] <- length(unique(result.coloc$results$snp))
  }
}
log_$sum.H3.H4 <- as.numeric(log_$sum.H3.H4)
log_ <- log_[order(log_$sum.H3.H4 , decreasing = T),]
log_$Coloc.Type[1] <- type_
log_$LD.Threshold[1] <- ld.threshold
log_$Distance[1] <- distance_
log_$UniqCpG.QTL[1] <- length(unique(QTLs$CPG))
log_$UniqSNPs.QTL[1] <- length(QTLs$SNP)
log_$UniqSNPs.GWAS[1] <- length(GWAS$SNP)
log_$CommonSNPs.QTL.GWAS[1] <- length(QTLs$SNP[QTLs$SNP %in% GWAS$SNP])

write.csv(log_,file=paste0(out.pref,".dist.",distance_,".coloc.",type_,".csv"),row.names=F)
save(result,file = paste0(out.pref,".dist.",distance_,".coloc.",type_,".rdat"))
