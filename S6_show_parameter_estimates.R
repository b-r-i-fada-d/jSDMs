# Set the base directory using your favorite method
# setwd("...")

# INPUT AND OUTPUT OF THIS SCRIPT
#	INPUT. the Fitted model.

# OUTPUT: Parameter estimates illustrated (for highest RUN of S2) in the file
# "results/parameter_estimates.pdf", the text file "results/parameter_estimates.txt", 
# as well as given numerically in multiple csv files (one per parameter type) named 
# "results/parameter_estimates_[parameter_name].csv".

#### --- Set up --- ####

set.seed(13)

localDir = "."
modelDir = file.path(localDir, "model")
resultDir = file.path(localDir, "results")


# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)

support.level.beta = NULL #Default: 0.95
support.level.gamma = NULL #Default: 0.95
support.level.omega = NULL #Default: 0.9
var.part.order.explained = NULL #Default: in variance partitioning of explained variance, species are shown in the order they are in the model
var.part.order.raw = NULL #Default: in variance partitioning of raw variance, species are shown in the order they are in the model
show.sp.names.beta = NULL #Default: species names shown in beta plot if there are at most 30 species and no phylogeny
plotTree = NULL #Default: tree is plotted in Beta plot if the model includes it
omega.order = NULL #Default: species shown in the order they are in the model
show.sp.names.omega = NULL #Default: species names shown in beta plot if there are at most 30 species


##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
#use support levels to selected the level of statistical support shown
#support.level.beta = 0.8
#support.level.gamma = 0.8
#support.level.omega = 0.8

#use var.part.order.explained to select which order species are shown in the raw variance partitioning
#var.part.order.raw should be a list of length the number of model. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "decreasing" if you wish to order according to explanatory power
#var.part.order.explained = list()
#var.part.order.explained[[1]] = 0
#var.part.order.explained[[2]] = c(2,1)

#use var.part.order.raw to select which order species are shown in the explained variance partitioning
#same options apply as for var.part.order.explained
#var.part.order.raw = list()
#var.part.order.raw[[1]] = "decreasing"
#var.part.order.raw[[2]] = c(1,2)

#use show.sp.names.beta to choose to show / not show species names in betaPlot
#if given, show.sp.names.beta should be a vector with length equalling number of model
#show.sp.names.beta = c(TRUE,FALSE)

#use plotTree to choose to plot / not plot the tree in betaPlot
#if given, plotTree should be a vector with length equalling number of model
#plotTree = c(FALSE,FALSE)

#use omega.order to select which order species are shown in omega plots
#omega.order should be a list of length the number of model. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "AOE" if you wish to use the angular order of the eigenvectors.
#omega.order = list()
#omega.order[[1]] = "AOE"
#omega.order[[2]] = c(2,1)
#Default: species shown in the order they are in the model
#show.sp.names.omega = c(TRUE,FALSE) #Default: species names shown in beta plot if there are at most 30 species
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################


if(is.null(support.level.beta)) support.level.beta = 0.95
if(is.null(support.level.gamma)) support.level.gamma =  0.95
if(is.null(support.level.omega)) support.level.omega =  0.9

library("Hmsc")
library("colorspace")
library("corrplot")
library("writexl")

samples_list = 500
thin_list = 1
nst = length(thin_list)
nChains = 4

text.file = file.path(resultDir,"/parameter_estimates.txt")
cat(c("This file contains additional information regarding parameter estimates.","\n","\n",sep=""),file=text.file)

filename = "model/m_FULL_500.Rdata" # EDIT HERE ###############
load(filename)
model <- m_FULL_500
rm(m_FULL_500)

# Append filename to text file
cat(c("\n", filename, "\n", "\n"), file = text.file, sep = "", append = TRUE)

# Initialize variables
nm = 1
var.part.order.explained = list(0)
var.part.order.raw = list(0)
omega.order = list(0)

# Extract model names
modelnames = "space"

# Open a PDF for plots
pdf(file = file.path(resultDir, "FULL_parameter_estimates.pdf"))


m = model

# Extract covariates
if (m$XFormula == "~.") {
  covariates = colnames(m$X)[-1]
} else {
  covariates = attr(terms(m$XFormula), "term.labels")
}

# Compute predictions and variance partitioning
preds = computePredictedValues(m)
VP = computeVariancePartitioning(m)
vals = VP$vals

# Evaluate model fit
MF = evaluateModelFit(hM = m, predY = preds)

# Store R2 values if available
R2 = MF$R2

# Save results to CSV
filename = file.path(resultDir, paste("FULL_parameter_estimates_VP_", ".csv", sep = ""))
write.csv(vals, file = filename)

# Save R2T Beta if available
if (!is.null(VP$R2T$Beta)) {
  filename = file.path(resultDir, paste("FULL_parameter_estimates_VP_R2T_Beta", modelnames[model_index], ".csv", sep = ""))
  write.csv(VP$R2T$Beta, file = filename)
}

# Save R2T Y if available
if (!is.null(VP$R2T$Y)) {
  filename = file.path(resultDir, paste("FULL_parameter_estimates_VP_R2T_Y", modelnames[model_index], ".csv", sep = ""))
  write.csv(VP$R2T$Y, file = filename)
}

# Plot variance partitioning
plotVariancePartitioning(hM = m, VP = VP, main = paste0("FULL_Proportion of explained variance, ", modelnames[model_index]))

# Close PDF output
dev.off()