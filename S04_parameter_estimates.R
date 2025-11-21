#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Post
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=2:00:00
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

# -----------------------------------------------------------------------------
# HMSC post-processing
# Assumptions:
#  - covariates>1 = YES  (model$nc > 1)
#  - traits>1     = NO   (model$nt == 1 or <=1) -> skip gamma plotting
#  - randomLevels>0 = YES (model$nr > 0)
# -----------------------------------------------------------------------------

# ---- outputs & model file ----
output_prefix <- "model"
if(is.null(arguments$output_dir) || arguments$output_dir == "") arguments$output_dir <- getwd()
dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)

parameter_file_path <- file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates.txt"))
cat("This file contains additional information regarding parameter estimates.\n\n", file = parameter_file_path)

# handle RDS naming convention (try _fitted.rds first)
model_rds_candidate <- gsub("\\.rds$", "_fitted.rds", arguments$model_rds)
model_rds_path <- if(file.exists(model_rds_candidate)) model_rds_candidate else arguments$model_rds
if(!file.exists(model_rds_path)) stop("Model file not found: ", model_rds_path)
message("Reading model from: ", model_rds_path)
model <- readRDS(model_rds_path)

# quick model diagnostics (write to parameter file)
cat("model dims and names:\n", file = parameter_file_path, append = TRUE)
cat("ns (species): ", model$ns, "\n", file = parameter_file_path, append = TRUE)
cat("nc (covariates): ", model$nc, "\n", file = parameter_file_path, append = TRUE)
cat("nt (traits): ", model$nt, "\n", file = parameter_file_path, append = TRUE)
cat("nr (random levels): ", model$nr, "\n\n", file = parameter_file_path, append = TRUE)

# open PDF for plots
plots_pdf <- file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates.pdf"))
pdf(file = plots_pdf)
message("Writing plots to: ", plots_pdf)

# -----------------------------------------------------------------------------
# VARIANCE PARTITIONING
# -----------------------------------------------------------------------------
message("Computing predicted values and variance partitioning...")
preds <- tryCatch(computePredictedValues(model), error = function(e) stop("computePredictedValues failed: ", e$message))
VP <- tryCatch(computeVariancePartitioning(model), error = function(e) stop("computeVariancePartitioning failed: ", e$message))
vals <- VP$vals
mycols <- grDevices::rainbow(nrow(vals))
MF <- tryCatch(evaluateModelFit(hM = model, predY = preds), error = function(e) { warning("evaluateModelFit failed: ", e$message); list() })

# attempt to pick an R2 metric (TjurR2, R2, SR2) if present
R2_candidates <- list(MF$TjurR2, MF$R2, MF$SR2)
R2 <- NULL
for(x in R2_candidates) if(!is.null(x)) { R2 <- x; break }
if(!is.null(R2)) vals <- rbind(vals, R2)
write.csv(vals, file = file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates_VP.csv")), row.names = FALSE)
plotVariancePartitioning(hM = model, VP = VP, main = "Proportion of explained variance", cex.main = 0.8, cols = mycols)

# If R2 available, plot raw variance partitioning scaled by R2
if(!is.null(R2)) {
  VPr <- VP
  # scale column k by R2[k]
  VPr$vals <- sweep(VPr$vals, 2, R2, FUN = "*")
  plotVariancePartitioning(hM = model, VP = VPr, main = "Proportion of raw variance", cex.main = 0.8, cols = mycols, ylim = c(0, 1))
}

# -----------------------------------------------------------------------------
# BETA: posterior estimates
# -----------------------------------------------------------------------------
message("Extracting Beta posterior estimates...")
postBeta <- getPostEstimate(model, parName = "Beta")
# defensive coercion to 2D (prevents collapse to vector)
postBeta$mean <- as.matrix(postBeta$mean)
postBeta$support <- as.matrix(postBeta$support)
postBeta$supportNeg <- as.matrix(postBeta$supportNeg)

# write Beta tables to xlsx (posterior mean, support positive, support negative)
me <- cbind(Species = model$spNames, as.data.frame(t(postBeta$mean)))
po <- cbind(Species = model$spNames, as.data.frame(t(postBeta$support)))
ne <- cbind(Species = model$spNames, as.data.frame(t(postBeta$supportNeg)))
beta_filename <- file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates_Beta.xlsx"))
tryCatch({ writexl::write_xlsx(list("Posterior mean" = me, "Pr(x>0)" = po, "Pr(x<0)" = ne), path = beta_filename) }, error = function(e) warning("Failed to write Beta xlsx: ", e$message))

# attempt to create a Beta plot
tryCatch({
  showNames <- is.null(model$phyloTree) && model$ns <= 100
  plotTreeHere <- !is.null(model$phyloTree)
  post <- postBeta
  plotBeta(model, post = post, supportLevel = 0.95, param = "Sign",
           plotTree = plotTreeHere,
           covNamesNumbers = c(TRUE, FALSE),
           spNamesNumbers = c(showNames, FALSE),
           cex = c(0.6, 0.6, 0.8))
  title(main = "BetaPlot", line = 2.5, cex.main = 0.8)
}, error = function(e) { warning("Beta plotting skipped/failed: ", e$message) })

# -----------------------------------------------------------------------------
# GAMMA: because traits>1 = NO (per user) we skip plotting Gamma but still
# extract posterior if present; coerce safely if it exists.
# -----------------------------------------------------------------------------
message("Handling Gamma (traits<=1 per user request: skipping plot)")
postGamma <- tryCatch(getPostEstimate(model, parName = "Gamma"), error = function(e) NULL)
if(!is.null(postGamma)) {
  postGamma$mean <- as.matrix(postGamma$mean)
  postGamma$support <- as.matrix(postGamma$support)
  postGamma$supportNeg <- as.matrix(postGamma$supportNeg)
  # write gamma to xlsx for inspection
  gamma_fn <- file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates_Gamma.xlsx"))
  tryCatch({
    me_g <- as.data.frame(t(postGamma$mean)); po_g <- as.data.frame(t(postGamma$support)); ne_g <- as.data.frame(t(postGamma$supportNeg))
    writexl::write_xlsx(list("Posterior mean" = me_g, "Pr(x>0)" = po_g, "Pr(x<0)" = ne_g), path = gamma_fn)
  }, error = function(e) warning("Failed to write Gamma xlsx: ", e$message))
}

# -----------------------------------------------------------------------------
# OMEGA / Associations: we assume randomLevels>0 per user instruction.
# process each random level WITHOUT for-loops (lapply) and coerce matrices
# -----------------------------------------------------------------------------
message("Computing associations (Omega) and writing outputs...")
OmegaCor <- computeAssociations(model)
# OmegaCor is a list; process each element with lapply
processOmega <- function(omegaListElem, rl_name) {
  # ensure mean and support are matrices
  mean_m <- as.matrix(omegaListElem$mean)
  support_m <- as.matrix(omegaListElem$support)
  
  # produce toPlot: cells that pass support threshold -> sign(mean)
  support_level <- 0.9
  toPlot <- ( (support_m > support_level) | (support_m < (1 - support_level)) ) * sign(mean_m)
  toPlot <- as.matrix(toPlot)
  
  # if toPlot has no dim or only one row/col, expand to square if possible
  if(is.null(dim(toPlot)) || nrow(toPlot) == 1 || ncol(toPlot) == 1) {
    # attempt to coerce to square using ns
    toPlot <- matrix(as.numeric(toPlot), nrow = model$ns, ncol = model$ns)
  }
  
  # remove species names if too many for plotting aesthetics
  showNamesOmega <- model$ns <= 100
  if(!showNamesOmega) {
    rownames(toPlot) <- rep("", nrow(toPlot))
    colnames(toPlot) <- rep("", ncol(toPlot))
  }
  
  # create corrplot (on current pdf device)
  main_title <- paste0("Associations: ", rl_name)
  tryCatch({
    corrplot(toPlot, method = "color",
             col = colorRampPalette(c("blue", "white", "red"))(20),
             mar = c(0, 0, 1, 0),
             main = main_title, cex.main = 0.8)
  }, error = function(e) warning("corrplot failed for ", rl_name, ": ", e$message))
  
  # write matrices to xlsx files
  me <- cbind(Species = c("", model$spNames), as.data.frame(mean_m)); colnames(me)[1] <- ""
  po <- cbind(Species = c("", model$spNames), as.data.frame(support_m)); colnames(po)[1] <- ""
  ne <- cbind(Species = c("", model$spNames), as.data.frame(1 - support_m)); colnames(ne)[1] <- ""
  out_fn <- file.path(arguments$output_dir, paste0(output_prefix, "_parameter_estimates_Omega_", rl_name, ".xlsx"))
  tryCatch({ writexl::write_xlsx(list("Posterior mean" = me, "Pr(x>0)" = po, "Pr(x<0)" = ne), path = out_fn) },
           error = function(e) warning("Failed to write Omega xlsx for ", rl_name, ": ", e$message))
  
  # also return a small summary for the text file
  summary_txt <- paste0("Omega (", rl_name, "): mean range [", round(min(mean_m, na.rm = TRUE), 3), ", ", round(max(mean_m, na.rm = TRUE), 3), "]\n")
  return(summary_txt)
}

omega_summaries <- mapply(function(o, nm) processOmega(o, nm), OmegaCor, names(OmegaCor), SIMPLIFY = TRUE)
# append omega summaries to parameter file
cat(paste(omega_summaries, collapse = "\n"), file = parameter_file_path, append = TRUE)

# -----------------------------------------------------------------------------
# finalization
# -----------------------------------------------------------------------------
message("Closing PDF device and finishing up...")
dev.off()
cat("Completed HMSC post-processing.\n", file = parameter_file_path, append = TRUE)
message("Done. Outputs in: ", normalizePath(arguments$output_dir))

# 
