# --- Set up --- ####
#	Input to script: fitted models
#	Output of script: MCMC convergence statistics for selected model parameters

library("Hmsc")
library("tidyverse")
library("ggplot2")
library("vioplot")
library("colorspace")

# set directories
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")

load(file.path(resultDir, "model_trial_fitted.RData"))

# make script reproducible
set.seed(13)

# --- Setting commonly adjusted parameters --- ####

showBeta = TRUE #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = TRUE #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = TRUE #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = 500 #Default: convergence of Omega shown for 50 randomly selected species pairs
showRho = TRUE #Default: showRho = FALSE, convergence not shown for rho-parameters
showAlpha = TRUE #Default: showAlpha = FALSE, convergence not shown for alpha-parameters


ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL
ma.rho = NULL
na.rho = NULL

# --- Define model parameters --- ####

samples = 100
thin = 2
nst = length(thin)
nChains = 4
Lst = 1

text.file = file.path(resultDir,"/MCMC_convergence.txt") # EDIT
cat("MCMC_Convergence_statistics\n\n",
    file <- text.file,
    sep = "") # EDIT

# --- Evaluate convergence --- ####

# filename = "models/m_SPACE_1000.RData"
# load(filename)
# models <- m_SPACE_1000 # EDIT
# rm(m_SPACE_1000)

# cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
nm <- 1
mpost <- convertToCodaObject(model, 
                             spNamesNumbers = c(T,F),
                             covNamesNumbers = c(T,F))
nr <- model$nr
cat(c("\n",names(model),"\n\n"),file=text.file,sep="",append=TRUE)
psrf <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
tmp <- summary(psrf)
cat("\nbeta\n\n", 
    file <- text.file,
    sep = "", 
    append=TRUE)

cat(tmp[,1], 
    file <- text.file, 
    sep = "\n", 
    append = TRUE)

ma.beta <- psrf[,1]
na.beta <- paste0(as.character(1), ",", as.character(thin)) # thin = 1; samples = 1000
ma.beta <- cbind(ma.beta,psrf[,1])
na.beta <- c(na.beta,paste0(as.character(1),",",as.character(samples))) # thin = 1; samples = 1000

psrf <- gelman.diag(mpost$Gamma, 
                    multivariate = FALSE)$psrf
tmp <- summary(psrf)
cat("\ngamma\n\n", 
    file <- text.file, 
    sep = "", 
    append=TRUE)
cat(tmp[,1], 
    file <- text.file, sep <- "\n", 
    append = TRUE)
ma.gamma <- psrf[,1]
na.gamma <- paste0(as.character(1), ",", as.character(thin)) # thin = 1; samples = 1000
ma.gamma <- cbind(ma.gamma, psrf[,1])
na.gamma <- c(na.gamma, paste0(as.character(1),",",as.character(samples))) # thin = 1; samples = 1000

# Write PSRF for Rho (if available)
if(showRho && !is.null(mpost$Rho)) {
  psrf <- gelman.diag(mpost$Rho, multivariate = FALSE)$psrf
  cat("\nrho\n\n",
      file <- text.file, 
      sep = "", 
      append = TRUE)
  cat(psrf[1], 
      file = text.file, 
      sep = "\n", 
      append = TRUE)
}

if(showOmega && nr > 0) {
  cat("\nomega\n\n", file = text.file, sep = "", append = TRUE)
  ma.omega <- NULL
  na.omega <- character(nr)
  for(k in 1:nr) {
    cat(c("\n", names(model$ranLevels)[k], "\n\n"), file = text.file, sep = "", append = TRUE)
    tmp <- mpost$Omega[[k]]
    z <- dim(tmp[[1]])[2]
    if(z > maxOmega) {
      sel <- sample(1:z, size = maxOmega)
      tmp <- lapply(tmp, function(mat) mat[, sel])
    }
    psrf <- gelman.diag(tmp, multivariate = FALSE)$psrf
    psrf_vals <- psrf[, 1]
    cat(psrf_vals, file = text.file, sep = "\n", append = TRUE)
    if(is.null(ma.omega)) {
      ma.omega <- matrix(psrf_vals, ncol = 1)
      na.omega[1] <- paste0(thin, ",", samples)
    } else {
      ma.omega <- cbind(ma.omega, psrf_vals)
      na.omega[k] <- paste0(thin, ",", samples)
    }
  }
  colnames(ma.omega) <- names(model$ranLevels)
  rownames(ma.omega) <- rownames(psrf)
}

if(showAlpha & nr>0){
        for(k in 1:nr){
          if(model$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n", 
                file <- text.file, 
                sep = "\n",
                append = TRUE)
            cat(c("\n",names(model$ranLevels)[k],"\n\n"), 
                file <- text.file, 
                sep = "",
                append = TRUE)
            psrf <- gelman.diag(mpost$Alpha[[k]], 
                                multivariate = FALSE)$psrf
            cat(psrf[,1],
                file <- text.file, 
                sep = "\n",
                append = TRUE)            
          }
        }
      }

  Lst = Lst + 1

pdf(file = file.path(resultDir,"/MCMC_convergence.pdf")) # EDIT
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta, main="psrf(beta)")
  legend("topright", legend="FULL", fill=rainbow_hcl(1))
  vioplot(ma.beta, ylim=c(0.9,1.1), main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma, main="psrf(gamma)")
  legend("topright", legend="Null_1")
  vioplot(ma.gamma, ylim=c(0.9,1.1), main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega, names = na.omega, main = "psrf(omega)")
  legend("topright", legend = "Null_2", fill = rainbow_hcl(1))
  vioplot(ma.omega, names = na.omega, ylim = c(0.9,1.1), main = "psrf(omega)")
}

dev.off()
