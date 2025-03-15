### MODELING OF MICROBES RELATED TO ANDROPOGON GERARDI: CONSTRUCT MODEL(S)
### Erica Newman & Adam B. Smith | 2024-07
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/S1_define_models.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/S1_define_models.r')
###
###	SCRIPT INPUT: Fit model(s)
### SCRIPT OUTPUT: MCMC convergence statistics for selected model parameters, illustrated (for all RUNs performed thus far in S3) in the file './results/MCMC_convergence.pdf', and the text file './results/MCMC_convergence.txt' (until MCMC convergence or computational limit is reached).

### setup
#########

rm(list = ls())
set.seed(1)

# workDir <- 'C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes'
workDir <- 'E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes'

setwd(workDir)

library(BayesLogit)
library(colorspace)
library(Hmsc)
library(omnibus)
library(vioplot)

### settings
############

set.seed(1)

showBeta <- NULL # default: showBeta <- TRUE, convergence shown for beta-parameters
showGamma <- NULL # default: showGamma <- FALSE, convergence not shown for gamma-parameters
showOmega <- NULL # default: showOmega <- FALSE, convergence not shown for Omega-parameters
maxOmega <- NULL # default: convergence of Omega shown for 50 randomly selected species pairs
showRho <- NULL # default: showRho <- FALSE, convergence not shown for rho-parameters
showAlpha <- NULL # default: showAlpha <- FALSE, convergence not shown for alpha-parameters

showBeta <- TRUE
showGamma <- TRUE
showOmega <- TRUE
maxOmega <- 100
showRho <- TRUE
showAlpha <- TRUE

dirCreate('./results')

if (is.null(showBeta)) showBeta = TRUE
if (is.null(showGamma)) showGamma = FALSE
if (is.null(showOmega)) showOmega = FALSE
if (is.null(maxOmega)) maxOmega = 50
if (is.null(showRho)) showRho = FALSE
if (is.null(showAlpha)) showAlpha = FALSE

convergence_file <- './results/MCMC_convergence.txt'

ma.beta <- NULL
na.beta <- NULL
ma.gamma <- NULL
na.gamma <- NULL
ma.omega<- NULL
na.omega <- NULL
ma.alpha <- NULL
na.alpha <- NULL  
ma.rho <- NULL
na.rho <- NULL

fit <- readRDS('./models/fit_model.rds')

mpost <- convertToCodaObject(fit, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
nr <- fit$nr

if (showBeta) {

	psrf <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
	tmp <- summary(psrf)
	
	say('Convergence:')
	print(tmp)
	
	if (ma.beta) {
		ma.beta = psrf[ , 1]
	} else {
		ma.beta <- cbind(ma.beta, psrf[ , 1])
		# na.beta <- c(na.beta, paste0(as.character(thin),',',as.character(samples)))
	}

}

if (showGamma) {

	psrf <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
	tmp <- summary(psrf)
	
	say('Convergence:')
	print(tmp)

	if(is.null(ma.gamma)){
		ma.gamma <- psrf[,1]
		na.gamma <- paste0(as.character(thin),', ',as.character(samples))
	} else {
		ma.gamma <- cbind(ma.gamma,psrf[,1])
		na.gamma <- c(na.gamma,paste0(as.character(thin),',',as.character(samples)))
	}


}

      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat('\ngamma\n\n',file=convergence_file,sep='',append=TRUE)
        cat(tmp[,1],file=convergence_file,sep='\n',append=TRUE)
        if(is.null(ma.gamma)){
          ma.gamma = psrf[,1]
          na.gamma = paste0(as.character(thin),',',as.character(samples))
        } else {
          ma.gamma = cbind(ma.gamma,psrf[,1])
          if(j==1){
            na.gamma = c(na.gamma,paste0(as.character(thin),',',as.character(samples)))
          } else {
            na.gamma = c(na.gamma,'')
          }
        }
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat('\nrho\n\n',file=convergence_file,sep='',append=TRUE)
        cat(psrf[1],file=convergence_file,sep='\n',append=TRUE)
      }
      if(showOmega & nr>0){
        cat('\nomega\n\n',file=convergence_file,sep='',append=TRUE)
        for(k in 1:nr){
          cat(c('\n',names(models[[j]]$ranLevels)[k],'\n\n'),file=convergence_file,sep='',append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=convergence_file,sep='\n',append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(thin),',',as.character(samples))
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(thin),',',as.character(samples)))
            } else {
              na.omega = c(na.omega,'')
            }
          }
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models[[j]]$ranLevels[[k]]$sDim>0){
            cat('\nalpha\n\n',file=convergence_file,sep='\n',append=TRUE)
            cat(c('\n',names(models[[j]]$ranLevels)[k],'\n\n'),file=convergence_file,sep='',append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=convergence_file,sep='\n',append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

pdf(file= file.path(resultDir,'/MCMC_convergence.pdf'))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main='psrf(beta)')
  legend('topright',legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main='psrf(beta)')
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main='psrf(gamma)')
  legend('topright',legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main='psrf(gamma)')
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main='psrf(omega)')
  legend('topright',legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main='psrf(omega)')
}
dev.off()

