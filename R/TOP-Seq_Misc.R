getSamples <- function(pathYAML = "../config_Samples.yml") {
    yaml::yaml.load_file(pathYAML)
}

attachInformation <- function(d, dInfo, 
    createDir = TRUE,
    pathDir = "../output/TOP-Seq_Analysis/data/reads/") {

    d$ID <- getID(d)
    d$pathFastq <- getPathFastq(d, dInfo$pathReads)

    if (createDir) {
        dir.create(paste0(pathDir, d$ID), recursive = TRUE, showWarnings = FALSE)
    }

    namesInfo <- names(dInfo)
    namesSample <- names(d)
    for(i in 1:length(dInfo)) {
        foo <- namesInfo[i]
        if (!foo %in% namesSample) {
            d[foo] <- dInfo[i]
        }
    }
    return(d)
}

getID <- function(d = sampleInfo) {
    if (any(names(d) == "ID_3")) {
        foo <- d$ID_3
    } else {
        foo <- grep("ID_[1-2]", names(d))
        foo <- paste(unlist(d[foo]), collapse = "@")
    }
    return(foo)
}

getPathFastq <- function(d = sampleInfo, pathDir) {
    list.files(file.path(pathDir, d$ID_1), 
               pattern = paste0(d$ID_2, ".*fastq"), 
               full.names = TRUE)
}
