fastqStat <- function(fileFastq, ID, dirOut = "../output/TOP-Seq_Analysis/results/fastqc/") {
    require(readFastQC)
    require(stringr)
    require(tidyverse)

    pathOutput <- paste0(dirOut, ID, "/")
    dir.create(pathOutput, recursive = TRUE, showWarnings = FALSE)
    fileResult <- paste0(pathOutput, 
                         gsub("\\.fastq", "", basename(fileFastq)), ".RDS")
    fileOutput <- paste0(pathOutput, 
                         gsub("\\.fastq$", "_fastqc", basename(fileFastq)))
    inFQ <- paste0(fileOutput, "/fastqc_data.txt")

    if (!file.exists(fileResult)) {
        system(paste("fastqc -o", pathOutput, "--nogroup --quiet", fileFastq))
    }
    unzip(paste0(fileOutput, ".zip"), exdir = pathOutput)

    res <- list()
    res[["nReads"]] <- subset(read_fastqc(inFQ), Measure == "Total Sequences")
    res[["BQ"]] <- read_fastqc(inFQ, module = "Per base sequence quality")
    res[["BC"]] <- read_fastqc(inFQ, module = "Per base sequence content")
    res[["SL"]] <- read_fastqc(inFQ, module = "Sequence Length Distribution")
    saveRDS(res, fileResult)
    unlink(fileOutput, TRUE)
    return(NULL)
}
runFASTQC <- function(sampleInfo) {
    files <- list.files(dirname(sampleInfo$pathFastq), 
                        pattern = "fastq$", full.names = TRUE)
    # Add original bam/fastq
    files <- c(getPathFastq(sampleInfo, dInfo$pathReads), files)

    foreach(i = files, .combine = rbind) %dopar% {
        fastqStat(i, sampleInfo$ID)
    }
    return(sampleInfo)
}
getReadNumber <- function(
    sampleInfo, 
    pathFastq = "../output/TOP-Seq_Analysis/results/fastqc/",
    pathMap  = "../output/TOP-Seq_Analysis/results/reads/Mapping_") {

    library(data.table)

    files <- list.files(paste0(pathFastq, sampleInfo$ID), pattern = "RDS",
        full.names = TRUE)
    res <- foreach(i = files, .combine = rbind) %do% {
        data.table(N = readRDS(i)$nReads$Value, 
                   Step = basename(i))
    }
    res[, Step := gsub("_.*", "", Step)]
    res[grep("^IonXpress", Step), Step := "Original"]

    nMapped <- fread(paste0("grep 'mapped' ", 
        pathMap, sampleInfo$ID, 
        ".log | head -n 1 | cut -f 1 -d ' '"))
    nHQ <- fread(paste0("wc -l ../output/TOP-Seq_Analysis/data/reads/", 
           sampleInfo$ID, "/", sampleInfo$ID, ".bed"))$V1
    nWithDup <- nrow(readRDS(gsub("noDup", "withDup", sampleInfo$pathNoDuplicates)))
    nNoDup <- nrow(readRDS(sampleInfo$pathNoDuplicates))
    nFromCG <- sum(readRDS(sampleInfo$pathCoverage)[, 3])

    res2 <- data.table(N = c(nMapped, nHQ, nWithDup, nNoDup, nFromCG),
                     Step = c("Mapped", "MappingQuality", "WithDup", "NoDup", "FromCG"))
    res <- rbind(res, res2)
    res[, ID := sampleInfo$ID]
    res[, Step := factor(Step, levels = c("Original", "LongReads", "Adapter5", 
        "Adapter3", "TrimQuality", "Mapped", "MappingQuality", 
        "WithDup", "NoDup", "FromCG"))][]
    return(res)
}
getDistanceToC <- function(
    sampleInfo, 
    pathDistance = "../output/TOP-Seq_Analysis/results/reads/distance_") {
    foo <- readRDS(paste0(pathDistance, sampleInfo$ID, ".RDS"))
    res <- sapply(0:10, function(x) mean(foo <= x))
    res <- c(res, mean(foo > 10))
    data.table(ID = sampleInfo$ID, distance = res, step = c(0:10, ">10"))
}
getIdentifiedCG <- function(
    sampleInfo, coverage = NULL, covGroup = "low") {
    library(data.table)

    if (is.null(coverage)) {
        if (covGroup == "low") {
            coverage <- seq(1:10)
        } else {
            coverage <- c(1, 10, 50, 100, 200)
        }
    }
    foo <- readRDS(sampleInfo$pathCoverage)[, 3]
    data.table(idenCG = sapply(coverage, function(x) sum(foo >= x)),
               coverage, ID = sampleInfo$ID)
}
getStatistics <- function(sampleInfo, coverageGroup = "low") {
    runFASTQC(sampleInfo)
    numberOfReads <- getReadNumber(sampleInfo)
    distances <- getDistanceToC(sampleInfo)
    identCG <- getIdentifiedCG(sampleInfo, covGroup = coverageGroup)
    covCG <- readRDS(sampleInfo$pathCoverage)[, 3]
    saveRDS(list(numberOfReads, distances,
                 identCG, covCG),
            paste0("../output/TOP-Seq_Analysis/results/reads/statistics_", 
                   sampleInfo$ID, ".RDS"))

}

