dPCA_main_impl <- function(nGroupNum, nDatasetNum, nSampleNum, nPaired, nLociNum, groupId, datasetId, repNum, sampleName, lociChr, lociStart, lociEnd, vdata, nTransform, nColMeanCent, nColStand, nMColMeanCent, nMColStand, dSNRCut, nUsedPCAZ, nUseRB, dPeakFDRCut) {
    .Call(irene_dPCA_main_impl, nGroupNum, nDatasetNum, nSampleNum, nPaired, nLociNum, groupId, datasetId, repNum, sampleName, lociChr, lociStart, lociEnd, vdata, nTransform, nColMeanCent, nColStand, nMColMeanCent, nMColStand, dSNRCut, nUsedPCAZ, nUseRB, dPeakFDRCut)
}

