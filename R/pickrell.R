norm.pickrell <- function(counts, x, sizeFactors = NULL, verbose = FALSE){
    ## Normaliztion as described in Pickrell et al (Nature, 2010)
    ##   Notation as in the supplementary material, page 10
    ## x : GC content
    ## counts : a count matrix
    bin.gc <- cut(x, unique(c(-1, quantile(x, seq(0.005, 1, 0.005)))))
    s <- apply(counts, 2, function(xx) tapply(xx, bin.gc, sum))
    f1 <- pmax(s, 0.5) / rowSums(s)
    if(is.null(sizeFactors)) {
        f0 <- colSums(s) / sum(s)
    } else {
        f0 <- sizeFactors / sum(sizeFactors)
    }
    f <- log2(f1) - log2(f0)
    bin.meangc <- tapply(x, bin.gc, mean)
    offset <- apply(f, 2, function(xx) {
        if(verbose) cat(".")
        fit1 <- smooth.spline(bin.meangc, xx)
        offset <- predict(fit1, x)$y
        offset
    }) ## This is \hat{g}_{ij}
    list(yc = counts * 2^(-offset), offset = offset)
}
