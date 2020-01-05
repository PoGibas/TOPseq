mapping <- function(
    sampleInfo = NULL,
    threads    = 40,
    pathDirRes = "../output/TOP-Seq_Analysis/results/reads/",
    pathWork   = "./") {

    outBAMtmp <- paste0(pathWork, "tmpBAM_", sampleInfo$ID, ".bam")
    outBAM <- paste0(dirname(sampleInfo$pathFastq), "/", sampleInfo$ID, ".bam")
    outBED <- paste0(dirname(sampleInfo$pathFastq), "/", sampleInfo$ID, ".bed")

    cmdMapping <- paste(
        "bwa mem",
            "-t", threads, 
            "-v 3",
            sampleInfo$pathGenome, sampleInfo$pathFastq, "|",
        "samtools view -b -S - >", outBAMtmp)
    cmdFlagstat <- paste0(
        "samtools flagstat ", outBAMtmp, " > ", 
        pathDirRes, "Mapping_", sampleInfo$ID, ".log")
    cmdSort <- paste(
        "samtools sort -@", threads, outBAMtmp, "|",
        "samtools view",
            "-b -@", threads, 
            "-q", sampleInfo$mappingQuality, 
            "-F 4 -o", outBAM , "-")
    cmdIndex <- paste("samtools index", outBAM)
    cmdBED <- paste("bedtools bamtobed -i", outBAM, ">", outBED)

    cmd <- paste(cmdMapping, cmdFlagstat, cmdSort,
                 cmdIndex, cmdBED, sep = " && ")
    if (!file.exists(outBED)) {
        system(cmd)
        unlink(c(outBAMtmp, paste0(outBAM, ".bai")))
    }
    sampleInfo$cmd_mapping <- cmd
    return(sampleInfo)
}

rmDuplicates <- function(sampleInfo, duplicatesBy = "Adapter3") {
    require(data.table)
    require(magrittr)

    ID <- sampleInfo$ID
    pathBED <- paste0(dirname(sampleInfo$pathFastq), "/", ID, ".bed")
    pathFQ <- paste0(dirname(sampleInfo$pathFastq), 
                     "/", duplicatesBy, "_", ID, ".fastq")
    outTMP <- paste0("dStart_", ID, ".RDS")
    outDup <- gsub("\\.bed$", "_withDuplicates.RDS", pathBED)
    outNoDup <- gsub("\\.bed$", "_noDuplicates.RDS", pathBED)

    if (!file.exists(outDup)) {
        # Start positions
        dStart <- fread(pathBED) %>%
            .[, st := as.numeric(V2)] %>%
            .[V6 == "-", st := as.numeric(V3) - 1] %>%
            .[, .(chr = V1, start = st, strand = V6, read = paste0("@", V4))] %>%
            setkey(read)
        # Select reads with only one mapped position 
        foo <- dStart[, .N, read][N == 1, read]
        saveRDS(dStart[read %in% foo], outTMP, compress = FALSE)
        rm(dStart)

        # FASTQ length
        dMain <- paste("paste - - - - <", pathFQ, 
                       "| awk '{print ($1), length($2)}'") %>%
            fread() %>%
            setkey(V1) %>%
            .[V1 %in% foo] %>%
            setnames(c("read", "length")) %>%
            merge(readRDS(outTMP), "read") %>%
            setkey(chr, start, strand, length)
        saveRDS(dMain[, .(chr, start, end = start + 1, strand, read)], outDup, compress = FALSE)
        dMain[, .(read = head(read, 1)), 
                .(chr, start, end = start + 1, strand, length)] %>%
            .[, .(chr, start, end, strand, read)] %>%
            setkey(chr, start, end) %>%
            saveRDS(outNoDup, compress = FALSE)
        unlink(outTMP)
    }
    sampleInfo$duplicatesRemovedBy <- pathFQ
    sampleInfo$pathNoDuplicates <- outNoDup
    return(sampleInfo)
}

getCoverage <- function(
    sampleInfo, 
    pathOutput = "../output/TOP-Seq_Analysis/data/coverage_density/coverage_",
    pathResult = "../output/TOP-Seq_Analysis/results/reads/distance_",
    DISTANCE = NULL) {

    if (is.null(DISTANCE)) {
        DISTANCE <- sampleInfo$distanceToCG         
    }

    sampleInfo$pathCoverage <- paste0(pathOutput, DISTANCE, "_", sampleInfo$ID, ".RDS")

    if (!file.exists(sampleInfo$pathCoverage)) {

        dCG <- readRDS(sampleInfo$pathCG)
        dTS <- readRDS(sampleInfo$pathNoDuplicates)
        res <- bedtools(args = "-d -t \"all\"", input = list(a = dTS, b = dCG))
        rm(dTS)

        # Read only to one CG
        foo <- res[, .N, V5][N == 1, V5]
        res <- res[V5 %in% foo] %>%
            .[, .(V6, V7, V4, V9)] %>%
            .[V9 >= 0] %>%
            setnames(c("chr", "start", "strand", "distance"))

        saveRDS(res$distance, paste0(pathResult, sampleInfo$ID, ".RDS"), compress = FALSE)

        res[distance <= DISTANCE, .N, .(chr, start)] %>%
            merge(dCG, c("chr", "start"), all.y = TRUE) %>%
            .[is.na(N), N := 0] %>%
            .[, .(chr, start, N)] %>%
            setkey(chr, start) %>%
            setnames("N", paste0("cov_", sampleInfo$ID)) %>%
            saveRDS(sampleInfo$pathCoverage, compress = FALSE)
    }
    return(sampleInfo)
}
getCoverageStrand <- function(
    sampleInfo, 
    pathOutput = "../output/TOP-Seq_Analysis/data/coverage_density/strand_coverage_",
    DISTANCE = NULL) {

    if (is.null(DISTANCE)) {
        DISTANCE <- sampleInfo$distanceToCG         
    }

    sampleInfo$pathCoverage <- paste0(pathOutput, DISTANCE, "_", sampleInfo$ID, ".RDS")

    if (!file.exists(sampleInfo$pathCoverage)) {

        dCG <- readRDS(sampleInfo$pathCG)
        dTS <- readRDS(sampleInfo$pathNoDuplicates)
        res <- bedtools(args = "-d -t \"all\"", input = list(a = dTS, b = dCG))
        rm(dTS)

        # Read only to one CG
        foo <- res[, .N, V5][N == 1, V5]
        res[V5 %in% foo] %>%
            .[, .(V6, V7, V4, V9)] %>%
            .[V9 >= 0] %>%
            setnames(c("chr", "start", "strand", "distance")) %>%
            .[distance <= DISTANCE, .N, .(chr, start, strand)] %>%
            setkey(chr, start) %>%
            saveRDS(sampleInfo$pathCoverage, compress = FALSE)
    }
    return(sampleInfo)
}

getCoverage0C <- function(
    sampleInfo, 
    pathOutput = "../output/TOP-Seq_Analysis/data/coverage_density/coverage_0C_") {

    sampleInfo$pathCoverage0C <- paste0(pathOutput, sampleInfo$ID, ".RDS")

    if (!file.exists(sampleInfo$pathCoverage0C)) {
 
        dCG <- readRDS(sampleInfo$pathCG) %>%
            setkey(chr, start, end)
        readRDS(sampleInfo$pathNoDuplicates) %>%
            .[, .(chr, 
                  start = ifelse(strand == "+", start - 1, start),
                  end   = ifelse(strand == "+", end, end + 1))] %>%
            .[start > 0] %>%
            .[, .N, .(chr, start, end)] %>%
            setkey(chr, start, end) %>%
            merge(dCG, ., c("chr", "start", "end"), all.x = TRUE) %>%
                .[is.na(N), N := 0] %>%
                .[, .(chr, start, N)] %>%
                setkey(chr, start) %>%
                setnames("N", paste0("cov_0C_", sampleInfo$ID)) %>%
                saveRDS(sampleInfo$pathCoverage0C, compress = FALSE)

    }
    return(sampleInfo)
}

getDensity <- function(
    sampleInfo, 
    pathOutput = "../output/TOP-Seq_Analysis/data/coverage_density/density_",
    chromosomes = paste0("chr", c(1:22, "X", "Y"))) {

    sampleInfo$pathDensity <- paste0(pathOutput, sampleInfo$ID, ".RDS")

    if (!file.exists(sampleInfo$pathDensity)) {

        dCG <- readRDS(sampleInfo$pathCG) %>%
            .[chr %in% chromosomes] %>%
            setkey(chr, start, end)
        densityCG <- list()
        for(i in chromosomes) {
            densityCG[[i]] <- kdeCG(dCG[chr == i])
        }
        dTS <- readRDS(sampleInfo$pathCoverage) %>%
            setnames(3, "cov")
        foreach(i = chromosomes, .combine = rbind) %do% {
            data.table(dTS[chr == i, .(chr, start)],
                       kde(dTS[chr == i, .(start, cov)]) / densityCG[[i]])
        } %>%
            setnames("dens", paste0("dens_", sampleInfo$ID)) %>%
            saveRDS(sampleInfo$pathDensity, compress = FALSE)
    }
    return(sampleInfo)
}
kdeCG <- function(data, bw = 80, np = 2^21, krnl = "epanechnikov",
    expL = 5e3, expH = 5e3) {
    lo <- min(data$start) - expL
    hi <- max(data$start) + expH
    d <- density(
        data$start, bw, kernel = krnl, 
        from = lo, to = hi, n = np)
    d <- ksmooth(
        d$x, d$y, kernel = "normal", bandwidth = bw, x.points = data$start)
    stopifnot(all(data$start == d$x))
    data.table(densCG = d$y)
}
kde <- function(data, bw = 180, np = 2^21, krnl = "epanechnikov") {
    library(data.table)
    lo <- min(data$start) 
    hi <- max(data$start)
    w <- data$cov / sum(data$cov)
    d <- density(
        data$start, bw, kernel = krnl, weights = w, 
        from = lo, to = hi, n = np)
    d <- ksmooth(
        d$x, d$y, kernel = "normal", bandwidth = bw, x.points = data$start)
    stopifnot(all(data$start == d$x))
    data.table(dens = d$y)
}