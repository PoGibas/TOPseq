findMotif <- function(
  path_genome = "/scratch/store/annotations/HomoSapiens_hg19/genome/genome.fa",
  motif = c("CG", "GC", "GCG")
) {
  library(data.table)
  library(foreach)
  library(magrittr)
  library(ShortRead)
  fasta <- readDNAStringSet(path_genome)
  for (MOTIF in motif) {
    print(MOTIF)
    foreach(i = seq_along(fasta), .combine = rbind) %do% {
      matchPattern(MOTIF, mask(fasta[[i]], "N"), fixed = FALSE) %>%
        start(.) %>%
        data.table(start = .) %>%
        .[, .(chr = names(fasta)[i], start)]
    } %>%
      saveRDS(paste0("motif_", MOTIF, ".RDS"), compress = FALSE)
  }
}