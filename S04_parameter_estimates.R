#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Post
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
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

model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"


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
      file = parameter_file_path, append = TRUE)
  cat(model$spNames[spOrderExplained], 
      file = parameter_file_path, sep = "\n", append = TRUE)
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
        file = parameter_file_path, append = TRUE)
    cat(model$spNames[spOrderRaw], 
        file = parameter_file_path, sep = "\n", append = TRUE)
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
    
    cat("\nomega.order\n\n", file = parameter_file_path, append = TRUE)
    cat(model$spNames[plotOrder], file = parameter_file_path, sep = "\n", append = TRUE)
    
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

