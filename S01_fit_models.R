#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Prep
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:15:00
#SBATCH --partition=cpu
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")

parser <- build.s01.argparser()
arguments <- parser$parse_args()

library("Hmsc")
library("tidyverse")
library("ggplot2")
library("jsonify")

dir.create(arguments$output_dir)

if(is.null(arguments$pa_data)) {
  # make df with unique station numbers, species PA, & env params
  species <- read_csv(arguments$species_data)
  env <- read_csv(arguments$env_data)
  
  # match up
  env <- env %>% select(-Station)
  species <- species %>% mutate(lat = round(lat, 4), lon = round(lon, 4))
  
  df <- left_join(species, env, by = c("lat", "lon", "Year", "Month"))
  df <- df %>% drop_na()
  
  rm(species, env)
} else {
  df <- read_csv(arguments$pa_data)
}

df <- df %>% drop_na()

#### --- Prep the model --- ####

# --- Select species data

Y <- df %>%
  dplyr::select(10:143) # just species columns
Y = as.matrix(Y)

# --- Select environmental data

XData = data.frame(year = as.factor(df$year), month = as.factor(df$month),
                   o2 = df$o2, temp = df$temp,
                   ph = df$ph, depth = df$depth)

# --- Coordinate data

xy = as.matrix(cbind(df$lon, df$lat))
rownames(xy) = as.character(df$station)
colnames(xy) = c("x-coordinate", "y-coordinate")
xy <- unique(xy)


#### --- Design model --- ####

# --- Study design

studyDesign <- data.frame(station = as.factor(df$station), 
                          year = as.factor(df$year) 
                         )

# --- Random effect structure (hierarchical study design)

rL.station = HmscRandomLevel(sData = xy, sMethod = "NNGP")
rL.year = HmscRandomLevel(units = df$year) #25.09 add back year rL

ranLevels = list(station = rL.station,
                 year = rL.year)


# --- Regression model for environmental covariates

XFormula = ~ o2 + temp + ph + depth + month + year

#### --- Construct the model --- ####

if (arguments$model_type == "full"){
  model = Hmsc(Y = Y, XData = XData,
               XFormula = XFormula,
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = ranLevels
               )
  )
} else if (arguments$model_type == "environmental"){
  model = Hmsc(Y = Y, XData = XData, 
               XFormula = XFormula,
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(year = rL.year
                                ))
} else if (arguments$model_type == "spatial"){
  model = Hmsc(Y = Y, XData = XData, 
               XFormula = ~1, # this is the difference
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(station = rL.station,
                                year = rL.year))
} else {stop(paste("unrecognised model type:", arguments$model_type))}


# Run model
sampler = sampleMcmc(model,
                     # initPar = "fixed effects",
                     engine="HPC",
                     thin = arguments$n_thins, 
                     samples = arguments$n_samples, 
                     transient = arguments$n_transients,
                     nChains = arguments$n_chains
                     )

init_file_path = file.path(arguments$output_dir, arguments$gibbs_samples_prefix)
saveRDS(to_json(sampler), file = init_file_path)

saveRDS(model, file = file.path(arguments$output_dir, arguments$model_rds))
