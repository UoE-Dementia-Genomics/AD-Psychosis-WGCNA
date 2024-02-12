args<-commandArgs(trailingOnly = TRUE)

InDir<-args[1]
FilePrefix <- args[2]
pvalue_threshold <- as.numeric(args[3])
save_rdat <- ifelse(args[4]=="T",T,F)
save_csv <- ifelse(args[5]=="T",T,F)
#print(pvalue_threshold)
library(stringr)
setwd(InDir)

cis_file = paste0(FilePrefix,"_Sig_mQTL_cis_",pvalue_threshold,".csv")
trans_file = paste0(FilePrefix,"_Sig_mQTL_trans_",pvalue_threshold,".csv")
merge_file = paste0(FilePrefix,"_mQTL_pvalue_",pvalue_threshold,".RData")
for (i in 1:22) {
  if(file.exists(paste0(FilePrefix,"_matrixEQTL_chr",i,"/",paste0(FilePrefix,"_matrixEQTL_chr",i,".RData")))){
    print(paste("Reading QTL file Chr",i))
    load(paste0(FilePrefix,"_matrixEQTL_chr",i,"/",paste0(FilePrefix,"_matrixEQTL_chr",i,".RData")))
    assign(x = paste0("mQTL_chr",i),value = me)
    remove(me)
  }

}
rm(i)
variables=ls()[str_detect(ls(),pattern ="mQTL_chr" )]
mQTLs <- mget(variables)
if(save_rdat){
  print("Saving merged QTL file...")
  save(mQTLs,file = merge_file)
}

if(save_csv){
  print(paste("Extracting Cis and Trans QTL with Pvalue <",pvalue_threshold,"..."))
  sig_mqtls_cis <- vector(mode = "list",length = length(mQTLs))
  sig_mqtls_trans <- vector(mode = "list",length = length(mQTLs))
  for (i in 1:length(mQTLs)) {
    sig_mqtls_cis[[i]] <- mQTLs[[i]]$cis$eqtls[mQTLs[[i]]$cis$eqtls$pvalue	< pvalue_threshold,]
    sig_mqtls_trans[[i]] <- mQTLs[[i]]$trans$eqtls[mQTLs[[i]]$trans$eqtls$pvalue	< pvalue_threshold,]
  }
  
  sig_mqtls_cis_out <- do.call(rbind.data.frame,sig_mqtls_cis)
  sig_mqtls_trans_out <- do.call(rbind.data.frame,sig_mqtls_trans)
  
  names(sig_mqtls_cis_out)[2] <- "gene"
  names(sig_mqtls_trans_out)[2] <- "gene"
  
  print(paste("Saving Cis and Trans QTL with Pvalue <",pvalue_threshold,"..."))
  write.csv(sig_mqtls_cis_out,file = paste0(FilePrefix,"_Sig_cis_eQTL_",pvalue_threshold,".csv"),row.names = F)
  write.csv(sig_mqtls_trans_out,file =paste0(FilePrefix,"_Sig_trans_eQTL_",pvalue_threshold,".csv"),row.names = F)  
}

