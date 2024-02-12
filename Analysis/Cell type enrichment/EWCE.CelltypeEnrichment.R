args <- commandArgs(T)

dir_10x <- args[1]
pheno_file <- args[2]
TargetGenes_file <- args[3]
rep_ <- as.numeric(args[4])
n_cores <- as.numeric(args[5])
generate_ctd <- ifelse(args[6]=="T",T,F)
ctd_file <- args[7]
out_prefix <- args[8]

cat("Input Arguments:\n")
cat("\n")
cat(paste("        10X files directory:",dir_10x,"\n"))
cat(paste("        Phenotype data:",pheno_file,"\n"))
cat(paste("        Target gene list:",TargetGenes_file,"\n"))
cat(paste("        Number of replication:",rep_,"\n"))
cat(paste("        Number of cores:",n_cores,"\n"))
cat(paste("        Generate Cell Type Dataset?",ifelse(generate_ctd,"Yes","No"),"\n"))
cat(paste("        CTD file:",ctd_file,"\n"))
cat(paste("        Output files prefix:",out_prefix,"\n"))
cat("\n")
cat("Loadning libraries...\n")
cat("\n")
suppressMessages(library(EWCE))
suppressMessages(library(Seurat))
suppressMessages(library(SummarizedExperiment))

set.seed(123456)

if(generate_ctd){
  cat("Reading Single Cell data...\n")
  cat("\n")
  sce <- Seurat::Read10X(dir_10x, gene.column = 1)
  pheno <- read.csv(pheno_file , stringsAsFactors = F)
  
  cat("Removing bad HGNC symbols...\n")
  cat("\n")
  sce1 <- fix_bad_hgnc_symbols(sce)
  
  cat("Droping uninformative genes...\n")
  cat("\n")
  sce2 <- drop_uninformative_genes(sce1, drop_nonhuman_genes=T , input_species = "human" ,output_species = "human" , level2annot = pheno$Subcluster, no_cores = n_cores)
  
  cat("Generating cell type dataset...\n")
  cat("\n")
  annotLevels = list(level1class=pheno$broad.cell.type,
                     level2class=pheno$Subcluster)
  ctd_file <- generate_celltype_data(exp=sce2,
                         annotLevels=annotLevels,
                         groupName=basename(out_prefix),
                         savePath=dirname(out_prefix), no_cores = n_cores) 
}
cat(paste("Loading cell type dataset from",ctd_file,"...\n"))
cat("\n")
ctd <- EWCE::load_rdata(ctd_file)


cat("Reading target genes...\n")
cat("\n")
TargetGenes <- read.csv(TargetGenes_file , stringsAsFactors = F , header = F)[,1]

cat("Cell type enrichment for level 1...\n")
cat("\n")
results1 <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = TargetGenes,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = rep_,
  annotLevel = 1, no_cores = n_cores)

cat("Cell type enrichment for level 2...\n")
cat("\n")
results2 <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = TargetGenes,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = rep_,
  annotLevel = 2, no_cores = n_cores)


cat("Saving results...\n")
cat("\n")
save(results1,results2,TargetGenes,ctd, file = paste0(out_prefix,".EWCE.rdat"))

cat("All done!\n")
