#' Returns upper triangular part of a matrix
#'
#' Used in sampling from `mnorm()`. Taken from Nimble manual section "5.2.4.1.2 LBJ distribution for correlation matrices" at https://r-nimble.org/manual/cha-writing-models.html#writing-models.
uppertri_mult_diag <- nimbleFunction(
    run = function(mat = double(2), vec = double(1)) {
        returnType(double(2))
        p <- length(vec)
        out <- matrix(nrow = p, ncol = p, init = FALSE)
        for(i in 1:p)
            out[ , i] <- mat[ , i] * vec[i]
        return(out)
})
