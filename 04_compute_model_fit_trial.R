# --- Set up --- ####

library("Hmsc")
library("tidyverse")

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
results.directory = file.path(localDir, "results")

set.seed(13)

# --- SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE --- ####

load("models/models_thin_1_samples_500_chains_4.Rdata")
model <- m_SPACE_500
rm(m_SPACE_500)

# thin = 1
# samples = 1000
# nChains = 4
# nfolds = NULL # Default: two-fold cross-validation
# cv.level = "year"
nParallel = 1 # Default: nParallel = nChains
 

preds = computePredictedValues(model)
MF = evaluateModelFit(hM = model, predY = preds)
# partition = createPartition(model, nfolds = nfolds, column = cv.level)
# predsCV = computePredictedValues(model,
#                                  partition = partition,
#                                  nParallel = nParallel
#                                  )
# MFCV = evaluateModelFit(hM = model, predY = predsCV)
WAIC = computeWAIC(model)
  
save(MF, 
     # MFCV, 
     WAIC, file = "results/SPACE_500.Rdata")
