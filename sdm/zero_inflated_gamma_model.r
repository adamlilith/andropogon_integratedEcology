p <- 0.2
theta <- 1
N <- 1000
x <- rexp(N, theta) * rbinom(N, 1, 1 - p)

dziexp <- nimbleFunction(run=function(x=double(1),p=double(0),theta=double(0),N=double(0),log=integer(0)){
  returnType(double(0))
  x0 <- x == 0
  return(sum(log(p * x0 + (1 - p) * dgamma(x, theta)*(!x0))))
})

rziexp <- nimbleFunction(run=function(n=integer(0),p=double(0),theta=double(0),N=double(0)){
  returnType(double(1))
  return(rexp(N,1)*rbinom(N,1,1-p))
})

# registerDistributions(list(
  # dziexp = list(
    # BUGSdist = 'dziexp(p,theta,N)',
    # types = c('value=double(1)','p=double(0)','theta=double(0)','N=double(0)'),
    # mixedSizes = TRUE
  # )
# ))

ziCode <- nimbleCode({
  x[1:n]~dziexp(p,theta,n)
  p~dunif(0,1)
  theta~dunif(0,10)
})

ziConsts <- list(n=N)
ziData <- list(x=x)
ziInits <- list(p=0.5,theta=1)

ziModel <- nimbleModel(ziCode, ziConsts,ziData,ziInits)
ziCModel <- compileNimble(ziModel,resetFunctions = TRUE)

zimcmcConf <- configureMCMC(ziModel,monitors=c('theta','p'),useConjugacy = FALSE)
zimcmcConf$removeSamplers('theta')
zimcmcConf$addSampler('theta',type='slice')

zimcmc<-buildMCMC(zimcmcConf)
ziCmcmc<-compileNimble(zimcmc, project = ziModel,resetFunctions = TRUE)
ziSamples<-runMCMC(ziCmcmc,
                      nchains = 2, nburnin=500,niter = 1000,samplesAsCodaMCMC = TRUE,
                      summary = FALSE, WAIC = FALSE)
plot(ziSamples)
