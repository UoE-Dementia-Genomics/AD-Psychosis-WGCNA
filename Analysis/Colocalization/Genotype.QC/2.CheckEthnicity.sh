
cd ${PROCESSDIR}

mkdir -p merge1KG_${FILEPREFIX}

# change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bim > merge1KG_${FILEPREFIX}/updateTo1KGFormat.txt
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --update-name merge1KG_${FILEPREFIX}/updateTo1KGFormat.txt --make-bed --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_QCd_1kgIDs

# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
plink --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_QCd_1kgIDs --bmerge ${RefDir}/1000g_GRCh37_HG19_Referance1.bed ${RefDir}/1000g_GRCh37_HG19_Referance1.bim ${RefDir}/1000g_GRCh37_HG19_Referance1.fam --maf 0.1 --geno 0.1 --make-bed --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G_test

## issue with variants at same position but different alleles - exclude these
plink --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_QCd_1kgIDs --exclude merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G_test-merge.missnp --make-bed --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_1kgIDs_forMerge

plink --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_1kgIDs_forMerge --bmerge ${RefDir}/1000g_GRCh37_HG19_Referance1.bed ${RefDir}/1000g_GRCh37_HG19_Referance1.bim ${RefDir}/1000g_GRCh37_HG19_Referance1.fam --maf 0.2 --geno 0.05 --make-bed --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G
# LD prune
plink --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G --indep 50 5 1.5 --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G.ld
plink --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G --extract merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G.ld.prune.in --make-bed --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G.ld.prune

rm merge1KG_${FILEPREFIX}/${FILEPREFIX}_1kgIDs_forMerge* 
rm merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G_test*

# use GCTA to calc PCs
gcta64 --bfile merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G.ld.prune --make-grm-bin --autosome --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G
gcta64 --grm merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G --pca --out merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G.pca

rm merge1KG_${FILEPREFIX}/${FILEPREFIX}_mergedw1000G*grm*

Rscript ${SCRIPTDIR}/Pre-Imputation/PlotEthnicity.R ${PROCESSDIR}  ${FILEPREFIX} ${RefDir} 

plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --remove merge1KG_${FILEPREFIX}/${FILEPREFIX}_EthnicityOutliers.txt --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd

echo "Check ethnithity is done"
