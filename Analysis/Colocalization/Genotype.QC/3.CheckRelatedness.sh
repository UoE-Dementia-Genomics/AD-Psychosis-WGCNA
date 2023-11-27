
cd ${PROCESSDIR}

## check for relatedness with other samples with KING
king -b QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bed --kinship --prefix QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd_king

## check for relatedness with other samples with plink
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --genome --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd_ibd

Rscript ${SCRIPTDIR}/Pre-Imputation/PlotRelatedness.r ${PROCESSDIR}/QCoutput_${FILEPREFIX}   ${FILEPREFIX}

plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --remove QCoutput_${FILEPREFIX}/${FILEPREFIX}_RelatednessOutliers.txt  --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd

echo "Relatedness was checked"


