##########################################################################
# This is an script to find genes near CpG locations using BioMart package
# The inputs are a list of CpGs, a manifest file to find the location of CpGs, a distance to define a searching window around each CpG, and the genome version (Default is "hg19")
##########################################################################


find.gene.list <- function(cpg.list , manifest, distance, genome){
  
  
  if(genome=="hg19"){
    ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",version = "GRCh37")
  }else
    if(genome=="hg38"){
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    }else 
      print("Please specify genome version: hg19 or hg38")
  
  filters <- c("chromosome_name","start","end")
  
  find.gene <- function(x){
    chr=as.character(x[2])
    start=as.character(as.numeric(x[3])-distance)
    end=as.character(as.numeric(x[4])+distance)
    values <- list(chromosome=chr,start=start,end=end)
    tryCatch({
      all.genes <- biomaRt::getBM(attributes=c("external_gene_name","ensembl_gene_id"), filters=filters, values=values, mart=ensembl)
      if(length(all.genes$external_gene_name)>0){
        genes <- paste(all.genes$external_gene_name,collapse =  ";")
        ens.id <- paste(all.genes$ensembl_gene_id,collapse =  ";")
        return(list(genes=genes,ens.id=ens.id))
      }
      else{
        return(list(genes="",ens.id=""))
      }
    },
    error = function(e) {return(list(genes="Error",ens.id="Error"))})
  }
  
  cpg <- data.frame(ID=cpg.list)
  cpg$chr = NA
  cpg$start = NA
  cpg$end = NA
  cpg$genes = NA
  cpg$ensembl = NA
  
  index = match(cpg$ID , manifest$IlmnID)
  cpg$chr = manifest$CHR[index]
  cpg$start = manifest$MAPINFO[index]
  cpg$end = manifest$MAPINFO[index]
  
  for (i in 1:nrow(cpg)) {
    results = find.gene(cpg[i,])
    cpg$genes[i] <- results$genes
    cpg$ensembl[i] <- results$ens.id
  }
  
  return(cpg)
}
#####################################################################################
library(vroom)

cpg_file = "CpG.List.csv" #simple txt or csv file. Each line contains a cpg ID
manifest_file= "MethylationEPIC_v-1-0_B4.csv" #Manifest file to find CpG position and chromosome
distance = 1000   #distance (bp) around CpG location to find genes
genome = "hg19"
out_pref = ""
  
cpg <- read.csv(file = cpg_file , header = T,stringsAsFactors = F)
manifest = vroom(manifest_file,skip = 7)

result <- find.gene.list(cpg.list = cpg$X,manifest = manifest,distance = distance,genome = genome)
write.csv(result,file = paste0(out_pref,".dist.",distance,".genes.csv"),row.names=F)

