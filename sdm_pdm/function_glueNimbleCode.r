#' Combine outputs of nimbleCode()
#'
#' @param ... Two or more outputs from [nimbleCode()].
#' 
#' @returns Concatenated `nimbleCode()` objects.
#'
#' @examples
#' a <- nimbleCode({x1 <- 23; x2 <- 13 })
#' b <- nimbleCode({x1 + x2})
#' c <- glueCode(a, b)
#' c
#' eval(c) # executes a and b together
glueNimbleCode <- function(...) {
    
    x <- list(...)
    
    # extract expressions from each object
    all_exprs <- unlist(
        lapply(x, function(expr) {
            if (is.call(expr) && identical(expr[[1]], as.name('{'))) {
                as.list(expr[-1])
            } else {
                list(expr)
            }
        }),
        recursive = FALSE
    )
    
    as.call(c(list(as.name('{')), all_exprs))
}
