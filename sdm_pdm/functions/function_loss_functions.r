#' MAE and RMSE loss functions to use in nimble::runCrossValidate()
#'
#' From nimble::runCrossValidate()
rmse_fx <- function(simulatedDataValues, actualDataValues){
  
    n <- length(simulatedDataValues) # simulatedDataValues is a vector
    sse <- 0
    for(i in 1:n){
        sse <- sse + (simulatedDataValues[i] - actualDataValues[i])^2
    }
    mse <- sse / n
    rmse <- sqrt(mse)
    rmse

}

#' Mean absolute error (better for skewed distributions)
mae_fx <- function(simulatedDataValues, actualDataValues){

    n <- length(simulatedDataValues) # simulatedDataValues is a vector
    ae <- 0
    for(i in 1:n){
      ae <- ae + abs(simulatedDataValues[i] - actualDataValues[i])
    }
    mae <- ae / n
    mae

}

#' Mean absolute percent error (better for skewed distributions)
mean_abs_percent_error_fx <- function(simulatedDataValues, actualDataValues){

    n <- length(simulatedDataValues) # simulatedDataValues is a vector
    pe <- 0
    for (i in 1:n) {
        pe <- pe + abs(simulatedDataValues[i] - actualDataValues[i]) / actualDataValues[i]
    }
    mpe <- pe / n
    mpe

}


