##############################################################################
##############################################################################

### This file fits the jSDM models ###
### for chapter 2 of BO'C PhD ###


# INPUT. Fitted models.

#	OUTPUT. Model fits computed by the cross-validation, with fitting
# (which is part of cross-validation) done for multiple RUNs

##############################################################################
##############################################################################

#### --- Set up --- ####

# load packages 
library("Hmsc")
library("tidyverse")
library("ggplot2")
library("colorspace")
library("vioplot")

# set directories
# setwd("C:/Users/boconnor/OneDrive - Marine Institute/Documents/Rwd/jSDMs/Sensitive-species-jSDMs")
localDir = "."
modelDir = file.path(localDir, "models")

set.seed(13)

##############################################################################
##############################################################################

# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE 

nfolds = NULL #Default: two-fold cross-validation
nParallel = 1 #Default: nParallel = nChains


##############################################################################
##############################################################################

if(is.null(nfolds)) nfolds = 2

samples_list = 800
thin_list = 1
nst = length(thin_list)
nChains = 4

if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename.in = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
  filename.out = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                          "_samples_", as.character(samples),
                                          "_chains_",as.character(nChains),
                                          "_nfolds_", as.character(nfolds),
                                          ".Rdata",sep = ""))
  if(file.exists(filename.out)){
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    print("model fit had been computed already")
  } else {
    if(file.exists(filename.in)){
      print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
      print(date())
      load(file = filename.in) #models
      nm = length(models)
      
      MF = list()
      MFCV = list()
      WAIC = list()
      
      for(mi in 1:nm){
        print(paste0("model = ",names(models)[mi]))
        m = models[[mi]]
        preds = computePredictedValues(m)
        MF[[mi]] = evaluateModelFit(hM=m, predY=preds)
        partition = createPartition(m, nfolds = nfolds) #USE column = ...
        preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
        MFCV[[mi]] = evaluateModelFit(hM=m, predY=preds)
        WAIC[[mi]] = computeWAIC(m)
      }
      names(MF)=names(models)
      names(MFCV)=names(models)
      names(WAIC)=names(models)
      
      save(MF,MFCV,WAIC,file = filename.out)
    }
  }
  Lst = Lst + 1
}
