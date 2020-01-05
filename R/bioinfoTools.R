writeBed <- function(data, name, ext = ".bed") {
  options(scipen = 999)
  data.table::fwrite(
    data, paste0(name, ext),
    quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE,
    nThread = 30
  )
}

bedtools <- function(
    method = "closest",
    args   = "",
    input  = list()) {
    require(magrittr)
    require(data.table)

    inptString <- ""
    fileBed <- c()
    for(i in 1:length(input)) {
        n <- names(input)[i]
        input[[i]] %>%
            setDT() %>%
            setkey(chr, start, end) %>%
            writeBed(n)
        inptString <- inptString %>%
            paste0(" -", n, " ", n, ".bed")
        fileBed <- c(fileBed, paste0(n, ".bed"))
    }
    res <- paste("bedtools", method, args, inptString) %>%
        fread(cmd = ., header = FALSE, nThread = 30)
 
    unlink(fileBed)
    return(res)
}

JC <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}
JCint <- function(x, y) {
    N <- length(intersect(x, y))
    data.frame(JC = N / length(union(x, y)), N)
}

getSRA <- function(SRA, name, pathRes = "./", fastq = TRUE, rmSRA = FALSE) {
  fileSRA <- downloadSRA(SRA, name, pathRes)
  if (fastq) {
    convertSRAtoFastq(fileSRA)
  }
  if (fastq & rmSRA) {
    unlink(fileSRA)
  }
  fileRes <- paste0(pathRes, "/", name, ifelse(fastq, ".fastq", ".sra"))
  return(fileRes)
}
downloadSRA <- function(SRA, name, pathRes, URL = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/", method = "internal", quiet = FALSE) {
  URLfull <- file.path(
    URL, substring(SRA, 1, 3), substring(SRA, 1, 6), SRA, paste0(SRA, ".sra")
  )
  fileSRA <- paste0(pathRes, "/", ifelse(is.null(name), SRA, name), ".sra")
  download.file(URLfull, fileSRA, method, quiet)
  return(fileSRA)
}
convertSRAtoFastq <- function(fileSRA, tool = "fastq-dump") {
  stopifnot(!is.null(fileSRA))
  if (!file.exists(fileSRA)) {
    stop(cat(fileSRA, "doesn't exist"))
  } else {
    system(paste(tool, "--outdir", dirname(fileSRA), fileSRA))
  }
  return(NULL)
}
