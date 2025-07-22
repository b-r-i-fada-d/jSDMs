#!/usr/bin/env Rscript

#SBATCH --jobname=HMSC-HPC.Prep
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")

parser <- build.s02.argparser()
arguments <- parser$parse_args()

library("Hmsc")
library("tidyverse")
library("ggplot2")
library("vioplot")
library("colorspace")
library("jsonify")

samples_list = vector("list", arguments$n_chains) # empty list with slots
for(i in 0:(arguments$n_chains-1)){
  chain_file <- arguments$gibbs_samples_prefix
  chain_file <- file.path(paste0(chain_file, "_chain_", i, ".rds"))
  samples_list[[i+1]] = from_json(readRDS(file = chain_file)[[1]])[[1]]
}

model <- readRDS(file.path(arguments$output_dir, arguments$model_rds))

fitSepTF = importPosteriorFromHPC(model,
                                  samples_list, 
                                  arguments$n_samples, 
                                  arguments$n_thins, 
                                  arguments$n_transients)

plotVariancePartitioning(fitSepTF, 
                         computeVariancePartitioning(fitSepTF), 
                         args.legend=list(x="bottomright"))

#### --- EVALUATE CONVERGENCE --- ####

showBeta = TRUE
showGamma = TRUE 
showOmega = TRUE 
showRho = TRUE 
showAlpha = TRUE
ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
maxOmega = 500
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL
ma.rho = NULL
na.rho = NULL

output_prefix <- tools::file_path_sans_ext(basename(arguments$model_rds))
convergence_file <- paste0(output_prefix, "_convergence.txt")
convergence_file <- file.path(arguments$output_dir, convergence_file) 
cat("MCMC_Convergence_statistics", file = convergence_file, sep = "")

mpost <- convertToCodaObject(model, 
                             spNamesNumbers = c(T,F), 
                             covNamesNumbers = c(T,F))
nr <- model$nr
cat("\nmodel\n\n", file = convergence_file, sep = "", append = TRUE)

# --- BETA: species interactions --- #
psrf_beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
cat("\nbeta\n\n", file = convergence_file, sep = "", append = TRUE)
cat(summary(psrf_beta)[,1], file = convergence_file, sep = "\n", append = TRUE)
ma.beta = psrf_beta[,1]
na.beta <- paste(arguments$thins, arguments$samples, sep = ",")

# --- GAMMA --- #
psrf_gamma = gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
cat("\ngamma\n\n", file = convergence_file, sep = "", append = TRUE)
cat(summary(psrf_gamma)[,1], file = convergence_file, sep = "\n", append = TRUE)
ma.gamma <- psrf_gamma[,1]
na.gamma <- paste(arguments$thins, arguments$samples, sep = ",")

# ---- RHO ---- #
if(showRho && !is.null(mpost$Rho)) {
  psrf_rho <- gelman.diag(mpost$Rho, multivariate = FALSE)$psrf
  cat("\nrho\n\n", file = convergence_file, sep = "", append = TRUE)
  cat(psrf_rho[1], file = convergence_file, sep = "\n", append = TRUE)
}

# --- OMEGA: random effects --- #
if(showOmega && nr > 0) {
  cat("\nomega\n\n", file = convergence_file, sep = "", append = TRUE)
  for (k in seq_len(nr)) {
    cat(paste0("\n", names(model$ranLevels)[k], "\n\n"), file = convergence_file, sep = "", append = TRUE)
    tmp <- mpost$Omega[[k]]
    z <- dim(tmp[[1]])[2]
    if(z > maxOmega) {
      sel <- sample(z, maxOmega)
      tmp <- lapply(tmp, function(x) x[, sel])
    }
    psrf_omega <- gelman.diag(tmp, multivariate = FALSE)$psrf
    cat(summary(psrf_omega)[,1], file = convergence_file, sep = "\n", append = TRUE)
    ma.omega <- cbind(ma.omega, psrf_omega[,1])
    na.omega <- c(na.omega, "1,1000")
  }
}
    
# --- ALPHA --- #
if(showAlpha && nr > 0) {
  for (k in seq_len(nr)) {
    if(model$ranLevels[[k]]$sDim > 0) {
      cat("\nalpha\n\n", file = convergence_file, sep = "\n", append = TRUE)
      cat(paste0("\n", names(model$ranLevels)[k], "\n\n"), file = convergence_file, sep = "", append = TRUE)
      psrf_alpha <- gelman.diag(mpost$Alpha[[k]], multivariate = FALSE)$psrf
      cat(psrf_alpha[,1], file = convergence_file, sep = "\n", append = TRUE)
    }
  }
}

# --- PLOT --- #

convergence_plots_path <- paste0(output_prefix, "_convergence.pdf")
convergence_plots_path <- file.path(arguments$output_dir, convergence_plots_path) 

pdf(file = convergence_plots_path)

if(showBeta) {
  par(mfrow = c(2,1))
  vioplot(ma.beta, col = rainbow_hcl(nm), names = na.beta, ylim = c(0, max(ma.beta)), main = "psrf(beta)")
  legend("topright", legend = "Beta", fill = rainbow_hcl(nm))
  vioplot(ma.beta, col = rainbow_hcl(nm), names = na.beta, ylim = c(0.9, 1.1), main = "psrf(beta)")
}

if(showGamma) {
  par(mfrow = c(2,1))
  vioplot(ma.gamma, col = rainbow_hcl(nm), names = na.gamma, ylim = c(0, max(ma.gamma)), main = "psrf(gamma)")
  legend("topright", legend = "Gamma", fill = rainbow_hcl(nm))
  vioplot(ma.gamma, col = rainbow_hcl(nm), names = na.gamma, ylim = c(0.9, 1.1), main = "psrf(gamma)")
}

if(showOmega && !is.null(ma.omega)) {
  par(mfrow = c(2,1))
  vioplot(ma.omega, col = rainbow_hcl(nm), names = na.omega, ylim = c(0, max(ma.omega)), main = "psrf(omega)")
  legend("topright", legend = "Omega", fill = rainbow_hcl(nm))
  vioplot(ma.omega, col = rainbow_hcl(nm), names = na.omega, ylim = c(0.9, 1.1), main = "psrf(omega)")
}

dev.off()

#### --- COMPUTE MODEL FIT --- ####

preds = computePredictedValues(model)
MF = evaluateModelFit(hM = model, predY = preds)
partition = createPartition(model,
                            nfolds = 5)

predsCV = computePredictedValues(model,
                                 partition = partition,
)
MFCV = evaluateModelFit(hM = model, predY = predsCV)
WAIC = computeWAIC(model)

fits_file <- paste0(output_prefix, "_fits.RData")
fits_file <- file.path(arguments$output_dir, fits_file)

save(MF, MFCV, WAIC, file = fits_file)


#### --- SHOW MODEL FIT --- ####

fits_plots_path <- paste0(output_prefix, "_fits.pdf")
fits_plots_path <- file.path(arguments$output_dir, fits_plots_path) 

pdf(file = fits_plots_path)

if(!is.null(MF$TjurR2)) {
  plot(MF$TjurR2, MFCV$TjurR2, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(model_label, ": Tjur R2\n",
                     "Mean(MF) = ", round(mean(MF$TjurR2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$TjurR2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}

if(!is.null(MF$R2)) {
  plot(MF$R2, MFCV$R2, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(model_label, ": R2\n",
                     "Mean(MF) = ", round(mean(MF$R2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$R2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}

if(!is.null(MF$AUC)) {
  plot(MF$AUC, MFCV$AUC, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(model_label, ": AUC\n",
                     "Mean(MF) = ", round(mean(MF$AUC, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$AUC, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0.5); abline(h = 0.5)
}

if(!is.null(MF$SR2)) {
  plot(MF$SR2, MFCV$SR2, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(model_label, ": SR2\n",
                     "Mean(MF) = ", round(mean(MF$SR2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$SR2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}
  
dev.off()


#### --- SHOW PARAMETER ESTIMATES --- ####

# Optional user-controlled parameters
support.level.beta <- NULL     # Default: 0.95
support.level.gamma <- NULL    # Default: 0.95
support.level.omega <- NULL    # Default: 0.9
var.part.order.explained <- NULL
var.part.order.raw <- NULL
show.sp.names.beta <- NULL
plotTree <- NULL
omega.order <- NULL
show.sp.names.omega <- NULL

### OUTPUTS ############################################################################

parameter_file_path <- paste0(output_prefix, "_parameter_estimates.txt")
parameter_file_path <- file.path(arguments$output_dir, parameter_file_path) 
cat("This file contains additional information regarding parameter estimates.\n\n",
    file = parameter_file_path, sep = "")

parameter_plots_path <- paste0(output_prefix, "_paramater_estimates.pdf")
parameter_plots_path <- file.path(arguments$output_dir, parameter_plots_path) 

pdf(file = parameter_plots_path)

# --- VARIANCE PARTITIONING --- #

if(model$XFormula == "~.") {
  covariates <- colnames(model$X)[-1]
} else {
  covariates <- attr(terms(model$XFormula), "term.labels")
}

if(model$nr + length(covariates) > 1 && model$ns > 1) {
  preds <- computePredictedValues(model)
  VP <- computeVariancePartitioning(model)
  vals <- VP$vals
  mycols <- rainbow(nrow(vals))
  MF <- evaluateModelFit(hM = model, predY = preds)
  
  R2 <- NULL
  if(!is.null(MF$TjurR2)) {
    vals <- rbind(vals, MF$TjurR2)
    R2 <- MF$TjurR2
  }
  if(!is.null(MF$R2)) {
    vals <- rbind(vals, MF$R2)
    R2 <- MF$R2
  }
  if(!is.null(MF$SR2)) {
    vals <- rbind(vals, MF$SR2)
    R2 <- MF$SR2
  }
  
  VP_file_path <- paste0(output_prefix, "_parameter_estimates_VP.csv")
  VP_file_path <- file.path(arguments$output_dir, VP_file_path)
  
  write.csv(vals, file = VP_file_path, row.names = FALSE)
    if(!is.null(VP$R2T$Beta))write.csv(VP$R2T$Beta, 
                                     file = file.path(arguments$output_dir, "parameter_estimates_VP_R2T_Beta.csv"))
  if(!is.null(VP$R2T$Y))write.csv(VP$R2T$Y, 
                                  file = file.path(arguments$output_dir, "parameter_estimates_VP_R2T_Y.csv"))
  
  # Explained variance order
  spOrderExplained <- if(all(var.part.order.explained == 0)) {
    seq_len(model$ns)
  } else if(all(var.part.order.explained == "decreasing")) {
    order(R2, decreasing = TRUE)
  } else {
    var.part.order.explained
  }
  
  VP$vals <- VP$vals[, spOrderExplained]
  cat("\nvar.part.order.explained\n\n", 
      file = text.file, append = TRUE)
  cat(model$spNames[spOrderExplained], 
      file = text.file, sep = "\n", append = TRUE)
  plotVariancePartitioning(hM = model, VP = VP, main = "Proportion of explained variance", cex.main = 0.8, cols = mycols)
  
  # Raw variance
  if(!is.null(R2)) {
    spOrderRaw <- if(all(var.part.order.raw == 0)) {
      seq_len(model$ns)
    } else if(all(var.part.order.raw == "decreasing")) {
      order(R2, decreasing = TRUE)
    } else {
      var.part.order.raw
    }
    
    VPr <- VP
    for (k in seq_len(model$ns)) {
      VPr$vals[, k] <- R2[k] * VPr$vals[, k]
    }
    VPr$vals <- VPr$vals[, spOrderRaw]
    cat("\nvar.part.order.raw\n\n", 
        file = text.file, append = TRUE)
    cat(model$spNames[spOrderRaw], 
        file = text.file, sep = "\n", append = TRUE)
    plotVariancePartitioning(hM = model, VP = VPr, main = "Proportion of raw variance", cex.main = 0.8, cols = mycols, ylim = c(0, 1))
  }
}

# --- BETA PLOTS --- #

if(model$nc > 1) {
  postBeta <- getPostEstimate(model, parName = "Beta")
  
  beta_filename <- paste0(output_prefix, "_parameter_estimates_Beta.xlsx")
  beta_filename <- file.path(arguments$output_dir, beta_filename)
  
  me <- cbind(Species = model$spNames, as.data.frame(t(postBeta$mean)))
  po <- cbind(Species = model$spNames, as.data.frame(t(postBeta$support)))
  ne <- cbind(Species = model$spNames, as.data.frame(t(postBeta$supportNeg)))
  
  writexl::write_xlsx(
    list("Posterior mean" = me, "Pr(x>0)" = po, "Pr(x<0)" = ne),
    path = beta_filename
  )
}

  
  showNames <- if(is.null(show.sp.names.beta))is.null(model$phyloTree) && model$ns <= 30 else show.sp.names.beta
  plotTreeHere <- !is.null(model$phyloTree) && (is.null(plotTree) || plotTree)
  
  plotBeta(model, post = postBeta, supportLevel = support.level.beta, param = "Sign",
           plotTree = plotTreeHere,
           covNamesNumbers = c(TRUE, FALSE),
           spNamesNumbers = c(showNames, FALSE),
           cex = c(0.6, 0.6, 0.8))
  mymain <- "BetaPlot"
  if(!is.null(model$phyloTree)) {
    rhovals <- unlist(poolMcmcChains(convertToCodaObject(model)$Rho))
    mymain <- paste0(mymain, ", E[rho] = ", round(mean(rhovals), 2), ", Pr[rho>0] = ", round(mean(rhovals > 0), 2))
  }
  title(main = mymain, line = 2.5, cex.main = 0.8)
}

# --- GAMMA PLOTS --- #

if(model$nt > 1 && model$nc > 1) {
  postGamma <- getPostEstimate(model, parName = "Gamma")
  plotGamma(model, post = postGamma, supportLevel = support.level.gamma, param = "Sign",
            covNamesNumbers = c(TRUE, FALSE),
            trNamesNumbers = c(model$nt < 21, FALSE),
            cex = c(0.6, 0.6, 0.8))
  title(main = "GammaPlot", line = 2.5, cex.main = 0.8)
}

# --- OMEGA PLOTS --- #

if(model$nr > 0 && model$ns > 1) {
  OmegaCor <- computeAssociations(model)
  for (r in seq_len(model$nr)) {
    toPlot <- ((OmegaCor[[r]]$support > support.level.omega) +
                 (OmegaCor[[r]]$support < (1 - support.level.omega)) > 0) *
      sign(OmegaCor[[r]]$mean)
    
    showNames <- if(is.null(show.sp.names.omega)) model$ns <= 30 else show.sp.names.omega
    if(!showNames) colnames(toPlot) <- rownames(toPlot) <- rep("", model$ns)
    
    plotOrder <- if(all(omega.order == 0)) {
      seq_len(model$ns)
    } else if(all(omega.order == "AOE")) {
      corrMatOrder(OmegaCor[[r]]$mean, order = "AOE")
    } else {
      omega.order
    }
    
    cat("\nomega.order\n\n", file = text.file, append = TRUE)
    cat(model$spNames[plotOrder], file = text.file, sep = "\n", append = TRUE)
    
    mymain <- paste0("Associations: ", names(model$ranLevels)[r])
    if(model$ranLevels[[r]]$sDim > 0) {
      alphavals <- unlist(poolMcmcChains(convertToCodaObject(model)$Alpha[[1]][, 1]))
      mymain <- paste0(mymain, ", E[alpha1] = ", round(mean(alphavals), 2), ", Pr[alpha1>0] = ", round(mean(alphavals > 0), 2))
    }
    
    corrplot(toPlot[plotOrder, plotOrder], method = "color",
             col = colorRampPalette(c("blue", "white", "red"))(3),
             mar = c(0, 0, 1, 0),
             main = mymain, cex.main = 0.8)
    
    me <- cbind("", as.data.frame(OmegaCor[[r]]$mean)); colnames(me)[1] <- ""
    po <- cbind("", as.data.frame(OmegaCor[[r]]$support)); colnames(po)[1] <- ""
    ne <- cbind("", 1 - as.data.frame(OmegaCor[[r]]$support)); colnames(ne)[1] <- ""
    vals <- list("Posterior mean" = me, "Pr(x>0)" = po, "Pr(x<0)" = ne)
    
    omega_filename <- paste0(output_prefix, "_parameter_estimates_Omega_", names(model$ranLevels)[r],".xlsx")
    omega_filename <- file.path(arguments$output_dir, omega_filename)
    writexl::write_xlsx(vals, path = omega_filename)
  }
}

dev.off()

