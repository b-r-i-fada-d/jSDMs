
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
#	INPUT. Model fits.
#	OUTPUT. Model fits illustrated (for highest RUN of S4) in the file "results/model_fit.pdf".

#### --- Set up --- ####

set.seed(13)

localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")

load("models/FULL_500_computed_fit.Rdata") # change this 

nfolds <- 2
samples <- 500
thin <- 1
nChains <- 4

nm <- length(MF)
modelnames <- names(MF)

pdf(file = "FULL_plots.pdf") # change this

if (!is.null(MF$TjurR2)) {
  plot(MF$TjurR2, main = "TjurR2", xlab = "Index", ylab = "Value")
}
if (!is.null(MF$AUC)) {
  plot(MF$AUC, main = "AUC", xlab = "Index", ylab = "Value")
}
if (!is.null(MF$RMSE)) {
  plot(MF$RMSE, main = "RMSE", xlab = "Index", ylab = "Value")
}

dev.off()
