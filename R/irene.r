dPCA <- function (meta, bed, data, minlen = 50, maxlen = 2e4, lambda = 0.25, verbose = TRUE, 
    nPaired = 0, nTransform = 0, nColMeanCent = 0, nColStand = 0, nMColMeanCent = 1, 
    nMColStand = 0, dSNRCut = 5, nUsedPCAZ = 0, nUseRB = 0, dPeakFDRCut = 0.5) 
{
    len <- abs(bed[, 3] - bed[, 2]) + 1
    groupId <- as.numeric(as.factor(meta$condition))
    datasetId <- as.numeric(as.factor(meta$factor))
    condId <- unique(groupId)
    nGroupNum <- length(condId)
    nDatasetNum <- length(unique(datasetId))
    nSampleNum <- nrow(meta)
    sampleName <- rep("", nSampleNum)
    repNum <- matrix(unlist(lapply(condId, function(i) {
        x <- datasetId[groupId == i]
        unlist(lapply(unique(x), function(j) {
            length(which((x == j)))
        }))
    })), nrow = nGroupNum, byrow = TRUE)
    data <- data[len > minlen, ] * 1e3/len
    if (min(data) <= 0) 
        data <- data - min(data) + 1e-5
    data <- normalize.quantiles(boxcox(data, lambda))
    if (verbose) 
        boxplot(data, xlab="datasets", ylab="boxcox(IP)",main="Normalized data")
    nLociNum <- nrow(bed)
    d <- dPCA_main_impl(nGroupNum, nDatasetNum, nSampleNum, 
            nPaired, nLociNum, groupId, datasetId, 
            repNum, sampleName, bed[, 1], bed[, 2], bed[, 3], 
            data, nTransform, nColMeanCent, nColStand, nMColMeanCent, 
            nMColStand, dSNRCut, nUsedPCAZ, nUseRB, dPeakFDRCut)
    d$proj <- matrix(unlist(lapply(1:4, function(i) d[[i]])), 
            nrow = 4, byrow = TRUE, dimnames = list(names(d)[1:4]))
    d$Dobs <- matrix(d$Dobs, ncol = nDatasetNum, byrow = TRUE)
    d$PC <- matrix(d$PC, ncol = nDatasetNum, byrow = TRUE)
    d
}

boxcox <- function(x,lambda=0.15) (x^lambda - 1)/lambda
