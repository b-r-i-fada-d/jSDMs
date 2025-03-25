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
model <- m_ENV_500
rm(m_ENV_500)

# thin = 1
# samples = 1000
# nChains = 4
nfolds = NULL # Default: two-fold cross-validation
cv.level = "station" 
nParallel = 1 # Default: nParallel = nChains
# if(is.null(nParallel)) nParallel = nChains
# if(is.null(nfolds)) nfolds = 2
# 
# 
# filename.in = file.path(model.directory, paste("models_thin_", as.character(thin), "_samples_", as.character(samples), "_chains_", as.character(nChains), ".Rdata", sep = ""))
# filename.out = model.directory/paste("MF_thin_", as.character(thin), "_samples_", as.character(samples), "_chains_", as.character(nChains), "_nfolds_", as.character(nfolds), "_cvlevel_", cv.level, ".Rdata", sep = "")
# 
# if (file.exists(filename.in)) {
#   print(paste0("thin = ", as.character(thin), "; samples = ", as.character(samples)))
#   print(date())
#   load(file = filename.in)
#   
#   model = m_ENV_1000

  preds = computePredictedValues(model)
  MF = evaluateModelFit(hM = model, predY = preds)
  partition = createPartition(model, nfolds = nfolds, column = cv.level)
  predsCV = computePredictedValues(model,
                                   partition = partition,
                                   nParallel = nParallel
                                   )
  MFCV = evaluateModelFit(hM = model, predY = predsCV)
  WAIC = computeWAIC(model)
  
  save(MF, 
       MFCV, 
       WAIC, file = "results/ENV_500_complete.Rdata")
# }
