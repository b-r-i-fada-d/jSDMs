# --- Set up --- ####

library("Hmsc")
library("tidyverse")
library("ggplot2")

# set directories
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)

# make script reproducible
set.seed(13)

#### --- Load data --- ####

# make df with unique station numbers, species PA, & env params

species <- read_csv("PA_species_thinned.csv")
env <- read_csv("stations_with_env__nosub_monthly.csv")

# match up
env <- env %>% select(-Station)
species <- species %>% mutate(lat = round(lat, 4), lon = round(lon, 4))

df <- left_join(species, env, by = c("lat", "lon", "Year", "Month"))

rm(species, env)

nas <- df %>% 
  filter(if_any(everything(), is.na)) # explore later

df <- df %>% drop_na()

##############################################################################
##############################################################################

#### --- Prep the model --- ####

# --- Select species data

Y <- df %>%
  dplyr::select(6:162) # just species columns
Y = as.matrix(Y)
hist(Y) # very zero inflated

# Explore species richness and prevalence

S = rowSums(Y)
P = colMeans(Y)
range(S)
range(P)
par(mfrow=c(1,2))
hist(S, xlab = "Species richness (S)")
hist(P, xlab = "Species prevalence (P)",xlim=c(0,1))


# --- Select environmental data

XData = data.frame(station = as.factor(df$Station), 
                   year = as.factor(df$Year), month = as.factor(df$Month),
                   o2 = df$o2, temp = df$temp,
                   ph = df$ph, depth = df$depth)

head(XData)
hist(XData$o2) # pretty normally distributed
hist(XData$temp)# pretty normally distributed
hist(XData$ph) # a bit left skewed
hist(XData$depth)

# --- Coordinate data

xy = as.matrix(cbind(df$lon, df$lat))
rownames(xy) = as.character(df$Station)
colnames(xy) = c("x-coordinate", "y-coordinate")
xy <- unique(xy)

head(xy)
plot(xy, asp=1) # show the map (NB., equal aspect ratio in the map)

##############################################################################
##############################################################################

#### --- Design model --- ####

# --- Study design

# Define the studyDesign as a dataframe 
# For example, if you have sampled the same locations over multiple years, you may define
# studyDesign = data.frame(sample = ..., year = ..., location = ...)
# studyDesign = data.frame(station=XData$station, depth=XData$depth)

studyDesign = data.frame(station = XData$station, year = XData$year)


# --- Random effect structure (hierarchical study design)

# For example, you may define year as an unstructured random effect
# rL.year = HmscRandomLevel(units = levels(studyDesign$year))

# For another example, you may define location as a spatial random effect
# rL.location = HmscRandomLevel(sData = locations.xy)
# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates

#rL = HmscRandomLevel(sData=xy)
rL.year = HmscRandomLevel(units = levels(studyDesign$year))
rL.station = HmscRandomLevel(units = levels(studyDesign$station))
# rL.depth = HmscRandomLevel(units = levels(studyDesign$depth))


# --- Regression model for environmental covariates

XFormula = ~ o2 + temp + ph + depth


# --- Construct the model

# Full model
m_FULL = Hmsc(Y = Y, XData = XData,  XFormula = XFormula,
         distr = "probit", # because PA
         studyDesign = studyDesign,
         ranLevels = list(station = rL.station,
                          year = rL.year))

# Spatial effect with no environmental variables
m_SPACE = Hmsc(Y = Y, XData = XData, 
               XFormula = ~1, # this is the difference
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(year = rL.year))

# Environmental covariates with no spatial random effect
m_ENV = Hmsc(Y = Y, XData = XData, 
             XFormula = XFormula,
             distr = "probit")

# --- Test that models fit without errors & save --- ####

m_FULL_test = sampleMcmc(m_FULL, samples = 2)
m_SPACE_test = sampleMcmc(m_SPACE, samples = 2)
m_ENV_test = sampleMcmc(m_ENV, samples = 2)

save(m_FULL_test, file = file.path(modelDir, "m_FULL_test_unfitted_model.RData"))
save(m_SPACE_test, file = file.path(modelDir, "m_SPACE_test_unfitted_model.RData"))
save(m_ENV_test, file = file.path(modelDir, "m_ENV_test_unfitted_model.RData"))

# all good!

# --- Fitting model --- ####

# We will store 100 posterior samples for each of two chains
# We note that for more "final" results, one might wish to have e.g. 1000 samples for each of four chains
# After fitting all models, we save the model object to a file
# Loading the fitted model then serves as the starting point for exploring the results
# The script runs over a loop where thin is first 1, then 10, then 100, and so on
# Thin is the thinning parameter of the MCMC chain.
# The transient (also called burn-in) is set to 50*thin
# When thin = 1, there will be 50 burn-in and 100 actual iterations. All actual iterations are stored.
# When thin = 10, there will be 500 burn-in and 1000 actual iterations. The actual iterations are thinned by 10, so 100 are stored.
# When thin = 100, there will be 5000 burn-in and 10000 actual iterations. The actual iterations are thinned by 100, so 100 are stored.
# A long MCMC chain is needed to achieve convergence
# Thinning is applied to avoid storing model objects of very large size
# Even if in the end thin = 1000 is required to achieve converge, We recommend to run the loop thin = 1, 10, 100, 1000
# This is for several reasons.
# First of all, it will not be known beforehand how much thinning is needed to achieve satisfactory convergence
# Second, thin = 1 will run very fast, whereas thin = 1000 will take very long (1000 times longer)
# After thin = 1 is completed, it is already possible to develop all the remaining scripts that explore the fitted model
# When exploring the fitted model, often one realizes changes that need to be made, even if the fitting has not converged
# Third, running the model fitting for thin = 1, 10, 100, 1000 does not take much longer than running it just for thin = 1000 (it takes ca. 12% longer)
# Thus, in summary, running the model fitting for thin = 1, 10, 100, 1000 typically saves a lot of time,
# as it allows one to proceed fast in writing (and revising) all the scripts that are needed from defining the model to producing the result tables and figures
# The idea is not to run the entire loop in one go, as that would take a lot of time. Just run thin = 1, and then move to develop the next scripts.
# You may then leave the remaining part of the loop (e.g. thin = 10, 100, 1000) to run e.g. overnight
# 
# nChains = 4
# nParallel = 2 # optional setting of nParallel
# samples = 1100
# thin = 1
# transient = 100
# 
# 
# # Run environmental model
# m_ENV_1000 = sampleMcmc(m_ENV, thin = thin, samples = samples, transient = transient,
#                         nChains = nChains, initPar = "fixed effects",
#                         nParallel = nParallel, verbose = verbose)
# 
# save(m_ENV_1000, file = file.path(modelDir, "m_ENV_1000_model.RData"))
# 
# # Run spatial model
# m_SPACE_1000 = sampleMcmc(m_SPACE, thin = thin, samples = samples, transient = transient,
#                         nChains = nChains, initPar = "fixed effects",
#                         nParallel = nParallel)
# 
# save(m_SPACE_1000, file = file.path(modelDir, "m_SPACE_1000_unfitted_model.RData"))
# 
# # Run full model
# m_FULL_1000 = sampleMcmc(m_FULL, thin = thin, samples = samples, transient = transient,
#                nChains = nChains, initPar = "fixed effects",
#                nParallel = nParallel)
# 
# #### --- Test & save models --- ####
# 
# models = list(m, m.ENV, m.ENV2, m.SPACE)
# names(models) = c("PA_model", "space_model", "env_model", "env2_model")
# save(models, file = file.path(modelDir, "unfitted_models.RData"))
# 
# # --- Testing that model fits without errors
# 
# for(i in 1:length(models)){
#   print(i)
#   sampleMcmc(models[[i]], samples = 2)
# }
# 
# save(m_FULL_1000, file = file.path(modelDir, "m_FULL_1000_unfitted_model.RData"))
# 
# # MCMC convergence can be difficult to achieve especially in those models that are not based on normal distribution.
# # For this reason, in the script above we initialize model with
# # initPar="fixed effects", with which option the MCMC chains are not started from locations randomized from the prior
# # but from a maximum likelihood solution to the fixed-effects part of the model
# 
# #################################################################################
# #################################################################################
# #################################################################################
# #################################################################################

#### --- Fit the models --- ####

models = list(m_FULL, m_ENV, m_SPACE)

# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE

nParallel = 1 # default: nParallel = nChains
# nParallel = 1 # set to 1 to disable parallel computing
# 
# load(file=file.path(modelDir,"unfitted_models.RData"))
nm = length(models)

samples_list = 1000
thin_list = 1
nChains = 4
transient = 100

if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    print("model had been fitted already")
  } else {
    print(date())
    for (mi in 1:nm) {
      print(paste0("model = ",names(models)[mi]))
      m = models[[mi]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     transient = transient,
                     nChains = nChains,
                     nParallel = nParallel) 
      models[[mi]] = m
    }
    save(models, file=filename)
  }
  Lst = Lst + 1
}

