catHeader <- function(text = "", level = 3) {
    cat(paste0("\n\n", 
               paste(rep("#", level), collapse = ""), 
               " ", text, "\n"))
}

fast_table <- function(x, y) {
  nx <- length(x)
  sx <- sum(x)
  sy <- sum(y)
  notx <- nx - sx
  tt <- sum(x & y)
  tf <- sx - tt
  ft <- sy - tt  
  ff <- notx - ft
  matrix(c(ff, tf, ft, tt), 2)
}