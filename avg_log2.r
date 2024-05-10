avg_log2 <- function(expr_array) {
    log2FC <- log(mean(expr_array), 2)
    return(log2FC)
}
