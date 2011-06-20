cqn.fixedlength <- function(counts, x, lengths, sizeFactors = NULL, subindex = NULL,
                            tau = 0.5, sqn = TRUE, verbose = FALSE){
    if(nrow(counts) != length(x) || nrow(counts) != length(lengths))
        stop("arguments 'counts' need to have the same number of rows as the length of arguments 'x' and 'lengths'")
    if(any(lengths <= 0))
        stop("argument 'lengths' need to be greater than zero")
    if(!is.null(sizeFactors) && length(sizeFactors) != ncol(counts))
        stop("argument 'sizeFactors' (when used) needs to have the same length as the number of columns of argument 'counts'")
    if(!is.null(subindex) && (min(subindex) <= 0 || max(subindex) > nrow(counts)))
        stop("argument 'subindex' (when used) needs to be indices into the number of rows of argument 'counts'")
    
    if(is.null(subindex))
        subindex <- which(rowMeans(counts) > 50)

    y <- sweep(sweep(log2(as.matrix(counts) + 1), 2, log2(sizeFactors/10^6)), 1, log2(lengths/1000))
    yfit <- y[subindex,, drop = FALSE]

    fixPredictor <- function(zz, varname = "B") {
        zz.fit <- zz[subindex]
        knots <- quantile(zz.fit, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) + c(+0.01, 0,0,0, -0.01)
        grid <- seq(from = min(zz.fit), to = max(zz.fit), length.out = 101)
        zz.out <- zz
        zz2 <- zz[-subindex]
        zz2[zz2 < min(zz.fit)] <- min(zz.fit)
        zz2[zz2 > max(zz.fit)] <- max(zz.fit)
        zz.out[-subindex] <- zz2
        list(fit = zz.fit, out = zz.out, knots = knots, grid = grid)
    }

    x1p <- fixPredictor(x)
    x1.fit <- x1p$fit
    x1.out <- x1p$out
    knots <- x1p$knots
    grid <- x1p$grid

    sp.fit <- ns(x1.fit, knots = knots)
    sp.predict <- predict(sp.fit, newx = x1.out)
    sp.func <- predict(sp.fit, newx = grid)

    if(verbose) cat("RQ fit ")
    regr <- lapply(1:ncol(yfit), FUN = function(ii) {
        if(verbose) cat(".")
        ytmp <- yfit[,ii]
        fit <- rq(ytmp ~ ., data = sp.fit, tau = tau)
        fitted <- predict(fit, newdata = sp.predict)
        func <- predict(fit, newdata = sp.func)
        list(fitted = fitted, func = func, coef = coef(fit))
    })
    if(verbose) cat("\n")
    
    fitted <- do.call(cbind, lapply(regr, function(xx) xx$fitted))
    func <- do.call(cbind, lapply(regr, function(xx) xx$func))
    ## This is the across sample median for the fitted values for the
    ## median value of the predictor amongst the subindex
    k <- order(x[subindex])[length(subindex) / 2]
    offset0 <- median(fitted[subindex[k],]) 

    residuals <- y - fitted
    if(sqn) {
        if(verbose) cat("SQN ")
        sqnFit <- SQN2(residuals, ctrl.id = subindex, min.q = 0.0001)
        residualsSQN <- sqnFit$yout
        if(verbose) cat(".\n")
        offset <- residualsSQN + offset0 - y
    } else {
        offset <- residuals + offset0 - y
    }
    
    rownames(offset) <- rownames(y)
    colnames(offset) <- colnames(y)
    colnames(func) <- colnames(y)
    
    out <- list(counts = counts, lengths = lengths, sizeFactors = sizeFactors, subindex = subindex, 
                y = y, x = x, offset = offset, offset0 = offset0,
                func1 = func, grid1 = grid, knots1 = knots)
    class(out) <- "cqn"
    out
}



cqn <- function(counts, x, lengths, sizeFactors = NULL, subindex = NULL,
                 tau = 0.5, sqn = TRUE, verbose = FALSE){
    if(nrow(counts) != length(x) || nrow(counts) != length(lengths))
        stop("arguments 'counts' need to have the same number of rows as the length of arguments 'x' and 'lengths'")
    if(any(lengths <= 0))
        stop("argument 'lengths' need to be greater than zero")
    if(!is.null(sizeFactors) && length(sizeFactors) != ncol(counts))
        stop("argument 'sizeFactors' (when used) needs to have the same length as the number of columns of argument 'counts'")
    if(!is.null(subindex) && (min(subindex) <= 0 || max(subindex) > nrow(counts)))
        stop("argument 'subindex' (when used) needs to be indices into the number of rows of argument 'counts'")
    
    if(is.null(subindex))
        subindex <- which(rowMeans(counts) > 50)
    y <- sweep(log2(as.matrix(counts) + 1), 2, log2(sizeFactors/10^6))
    yfit <- y[subindex,, drop = FALSE]

    fixPredictor <- function(zz, varname = "B") {
        zz.fit <- zz[subindex]
        knots <- quantile(zz.fit, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) + c(+0.01, 0,0,0, -0.01)
        grid <- seq(from = min(zz.fit), to = max(zz.fit), length.out = 101)
        zz.out <- zz
        zz2 <- zz[-subindex]
        zz2[zz2 < min(zz.fit)] <- min(zz.fit)
        zz2[zz2 > max(zz.fit)] <- max(zz.fit)
        zz.out[-subindex] <- zz2
        list(fit = zz.fit, out = zz.out, knots = knots, grid = grid)
    }

    x1p <- fixPredictor(x)
    x1.knots <- x1p$knots

    x2p <- fixPredictor(log2(lengths/1000))
    x2.knots <- x2p$knots

    df.fit <- data.frame(x1 = x1p$fit, x2 = x2p$fit)
    df.predict <- data.frame(x1 = x1p$out, x2 = x2p$out)
    df.func <- rbind(data.frame(x1 = x1p$grid,
                                x2 = median(x2p$grid)),
                     data.frame(x1 = median(x1p$grid),
                                x2 = x2p$grid))

    if(verbose) cat("RQ fit ")
    regr <- lapply(1:ncol(yfit), FUN = function(ii) {
        if(verbose) cat(".")
        ytmp <- yfit[,ii]
        fit <- rq(ytmp ~ ns(x1, knots = x1.knots) +
                  ns(x2, knots = x2.knots), data = df.fit, tau = tau)
        fitted <- predict(fit, newdata = df.predict)
        func <- predict(fit, newdata = df.func)
        list(fitted = fitted, func = func, coef = coef(fit))
    })
    if(verbose) cat("\n")
    
    fitted <- do.call(cbind, lapply(regr, function(xx) xx$fitted))
    func <- do.call(cbind, lapply(regr, function(xx) xx$func))
    func1 <- func[seq(from = 1, length.out = length(x1p$grid)),, drop = FALSE]
    func2 <- func[seq(from = length(x1p$grid) + 1, length.out = length(x2p$grid)),,drop = FALSE]
    
    ## This is the across sample median for the fitted values for the
    ## median value of the predictor amongst the subindex
    k <- order(x[subindex])[length(subindex) / 2]
    offset0 <- median(fitted[subindex[k],]) 

    residuals <- y - fitted
    if(sqn) {
        if(verbose) cat("SQN ")
        sqnFit <- SQN2(residuals, ctrl.id = subindex, min.q = 0.0001)
        residualsSQN <- sqnFit$yout
        if(verbose) cat(".\n")
        offset <- residualsSQN + offset0 - y
    } else {
        offset <- residuals + offset0 - y
    }
        
    rownames(offset) <- rownames(y)
    colnames(offset) <- colnames(y)
    colnames(func1) <- colnames(func2) <- colnames(y)

    out <- list(counts = counts, lengths = lengths, sizeFactors = sizeFactors, subindex = subindex, 
                y = y, x = x, offset = offset, offset0 = offset0, 
                func1 = func1, func2 = func2,
                grid1 = x1p$grid, grid2 = x2p$grid,
                knots1 = x1p$knots, knots2 = x2p$knots)
    class(out) <- "cqn"
    out
}

print.cqn <- function(x, ...) {
    cat("object of class 'cqn' with")
    cat(" ", nrow(x$counts), "regions\n")
    cat(" ", ncol(x$counts), "samples\n")
    cat(" originating from a fit with", ifelse(!is.null(x$func2), "smooth", "fixed (offset)"),
        "length\n")
    invisible(x)
}

cqnplot <- function(x, n = 1, col = "grey60", ylab = "QR fit", xlab = "", type = "l", lty = 1, ...) {
    if(class(x) != "cqn")
        stop("'x' needs to be of class 'cqn'")
    if(n == 1) {
        func <- x$func1
        grid <- x$grid1
        knots <- x$knots1
    }
    if(n == 2) {
        func <- x$func2
        grid <- x$grid2
        knots <- x$knots2
    }
    matplot(replicate(ncol(func), grid), func, ylab = ylab, xlab = xlab, type = type,
            col = col, lty = lty, ...)
    rug(knots, lwd = 2)
    invisible(x)
}

