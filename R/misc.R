catHeader <- function(text = "", level = 3) {
    cat(paste0("\n\n", 
               paste(rep("#", level), collapse = ""), 
               " ", text, "\n"))
}