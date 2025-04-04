# --- Set up --- ####

library("Hmsc")
library("tidyverse")

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
results.directory = file.path(localDir, "results")

set.seed(13)

# --- SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE --- ####

filename = "models/m_ENV_500.RData"
load(filename)
model <- m_ENV_500 # EDIT
rm(m_ENV_500)

thin = 1
samples = 500
nChains = 4
nfolds = 1 # Default: two-fold cross-validation
# cv.level = "year"
nParallel = 1 # Default: nParallel = nChains
 

preds = computePredictedValues(model)
MF = evaluateModelFit(hM = model, predY = preds)
partition = createPartition(model, nfolds = nfolds
                            # , column = cv.level
                            )

predsCV = computePredictedValues(model,
                                 partition = partition,
                                 nParallel = nParallel
                                 )
MFCV = evaluateModelFit(hM = model, predY = predsCV)
WAIC = computeWAIC(model)
  
save(MF, 
     MFCV, 
     WAIC, file = "results/model_fit_ENV_500.Rdata")
