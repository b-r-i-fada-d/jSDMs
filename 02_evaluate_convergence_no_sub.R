# --- Set up --- ####
#	Input to script: fitted models
#	Output of script: MCMC convergence statistics for selected model parameters

library("Hmsc")
library("tidyverse")
library("ggplot2")

# set directories
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)
resultDir = file.path(localDir, "results_june")
if(!dir.exists(resultDir)) dir.create(resultDir)

# make script reproducible
set.seed(13)

# --- Setting commonly adjusted parameters --- ####

showBeta = TRUE #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = TRUE #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = TRUE #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = 100 #Default: convergence of Omega shown for 50 randomly selected species pairs
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

samples_list = 1000
thin_list = 1
nst = length(thin_list)
nChains = 4
Lst = 1

text.file = file.path(resultDir,"/MCMC_convergence_no_sub_1000.txt") # EDIT
cat("MCMC_Convergence_statistics_no_sub_1000\n\n",file=text.file,sep="") # EDIT

# --- Evaluate convergence --- ####

filename = "models/nosub_models_thin_1_samples_1000_chains_4.Rdata"
load(filename)
models <- nosub_models_thin_1_samples_1000_chains_4 # EDIT
rm(nosub_models_thin_1_samples_1000_chains_4)

cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
nm = length(models)
mpost = convertToCodaObject(models, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
nr = models$nr
cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
tmp = summary(psrf)
cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
ma.beta = psrf[,1]
na.beta = paste0(as.character(1),",",as.character(1000)) # thin = 1; samples = 1000
ma.beta = cbind(ma.beta,psrf[,1])
na.beta = c(na.beta,paste0(as.character(1),",",as.character(1000))) # thin = 1; samples = 1000

psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
tmp = summary(psrf)
cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
ma.gamma = psrf[,1]
na.gamma = paste0(as.character(1),",",as.character(1000)) # thin = 1; samples = 1000
ma.gamma = cbind(ma.gamma,psrf[,1])
na.gamma = c(na.gamma,paste0(as.character(1),",",as.character(1000))) # thin = 1; samples = 1000

      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
        cat(psrf[1],file=text.file,sep="\n",append=TRUE)
      }

      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
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
          cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(1),",",as.character(1000)) # thin = 1; samples = 1000
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(1),",",as.character(1000))) # thin = 1; samples = 1000
            } else {
              na.omega = c(na.omega,"")
            }
          }
        }
      }
      
if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
            cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
          }
        }
      }

  Lst = Lst + 1

pdf(file= file.path(resultDir,"/MCMC_convergence_no_sub_1000.pdf")) # EDIT
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  # legend("topright", legend = "FULL", fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  # legend("topright",legend = "Null_1", fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  # legend("topright",legend = "Null_2", fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()
