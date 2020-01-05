filterLength <- function(
    sampleInfo, 
    pathDirReads = "../output/TOP-Seq_Analysis/data/reads/") {

    pathOutput <- paste0(pathDirReads, sampleInfo$ID, 
                         "/LongReads_", sampleInfo$ID, ".fastq")
    cmd <- paste("fastq_quality_trimmer", 
            "-t 1",
            "-l", sampleInfo$readLength, 
            "-i", sampleInfo$pathFastq,
            "-o", pathOutput,
            "-Q 33")
    if (!file.exists(pathOutput)) {
        system(cmd)
    }
    sampleInfo$pathFastq <- pathOutput
    sampleInfo$cmd_filterLength <- cmd
    return(sampleInfo)
}

cutadapt5 <- function(
    sampleInfo,
    pathDirReads = "../output/TOP-Seq_Analysis/data/reads/",
    pathDirRes   = "../output/TOP-Seq_Analysis/results/reads/") {

    file <- paste(c("Adapter5",
                    unlist(strsplit(basename(sampleInfo$pathFastq), "_"))[-1]),
                    collapse = "_")
    pathOutput <- paste0(pathDirReads, sampleInfo$ID, "/", file)
    cmd <- paste("cutadapt", 
            "-g", sampleInfo$adapter5Seq, 
            "-e", sampleInfo$adapter5Error,
            "-O", 10,
            "--trimmed-only",
            "-o", pathOutput, 
            sampleInfo$pathFastq, 
            "| tee", paste0(pathDirRes, gsub("fastq", "log", file)))
    if (!file.exists(pathOutput)) {
        system(cmd)
    }
    sampleInfo$pathFastq <- pathOutput
    sampleInfo$cmd_cutadapt5 <- cmd
    return(sampleInfo)
}

cutadapt3 <- function(
    sampleInfo,
    pathDirReads = "../output/TOP-Seq_Analysis/data/reads/",
    pathDirRes   = "../output/TOP-Seq_Analysis/results/reads/") {

    file <- paste(c("Adapter3",
                    unlist(strsplit(basename(sampleInfo$pathFastq), "_"))[-1]),
                    collapse = "_")
    pathOutput <- paste0(pathDirReads, sampleInfo$ID, "/", file)
    cmd <- paste("cutadapt", 
            "-a", sampleInfo$adapter3Seq, 
            "-e", sampleInfo$adapter3Error,
            "-O", 10,
            "-o", pathOutput, 
            sampleInfo$pathFastq, 
            "| tee", paste0(pathDirRes, gsub("fastq", "log", file)))
    if (!file.exists(pathOutput)) {
        system(cmd)
    }
    sampleInfo$pathFastq <- pathOutput
    sampleInfo$cmd_cutadapt3 <- cmd
    return(sampleInfo)
}

trimQuality <- function(
    sampleInfo,
    pathDirReads = "../output/TOP-Seq_Analysis/data/reads/") {

    file <- paste(c("TrimQuality",
                    unlist(strsplit(basename(sampleInfo$pathFastq), "_"))[-1]),
                    collapse = "_")
    pathOutput <- paste0(pathDirReads, sampleInfo$ID, "/", file)
    cmd <- paste("fastq_quality_trimmer", 
            "-t", sampleInfo$trimQuality,
            "-l", sampleInfo$trimLength, 
            "-i", sampleInfo$pathFastq, 
            "-o", pathOutput, 
            "-Q 33")
    if (!file.exists(pathOutput)) {
        system(cmd)
    }
    sampleInfo$pathFastq <- pathOutput
    sampleInfo$cmd_trimQuality <- cmd
    return(sampleInfo)
}