#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Post
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=05:00:00
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

#### --- COMPUTE MODEL FIT --- #####

model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"

preds = computePredictedValues(model)
MF = evaluateModelFit(hM = model, predY = preds)
partition = createPartition(model,
                            nfolds = 2)

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
       main = paste0(output_prefix, ": Tjur R2\n",
                     "Mean(MF) = ", round(mean(MF$TjurR2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$TjurR2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}

if(!is.null(MF$R2)) {
  plot(MF$R2, MFCV$R2, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(output_prefix, ": R2\n",
                     "Mean(MF) = ", round(mean(MF$R2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$R2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}

if(!is.null(MF$AUC)) {
  plot(MF$AUC, MFCV$AUC, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(output_prefix, ": AUC\n",
                     "Mean(MF) = ", round(mean(MF$AUC, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$AUC, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0.5); abline(h = 0.5)
}

if(!is.null(MF$SR2)) {
  plot(MF$SR2, MFCV$SR2, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Explanatory Power", ylab = "Predictive Power",
       main = paste0(output_prefix, ": SR2\n",
                     "Mean(MF) = ", round(mean(MF$SR2, na.rm = TRUE), 3),
                     ", Mean(MFCV) = ", round(mean(MFCV$SR2, na.rm = TRUE), 3)))
  abline(0, 1); abline(v = 0); abline(h = 0)
}

dev.off()
