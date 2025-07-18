#!/usr/bin/env Rscript

#SBATCH --jobname=HMSC-HPC.Prep
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
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

# make df with unique station numbers, species PA, & env params
species <- read_csv(arguments$species_data)
env <- read_csv(arguments$env_data)

# match up
env <- env %>% select(-Station)
species <- species %>% mutate(lat = round(lat, 4), lon = round(lon, 4))

df <- left_join(species, env, by = c("lat", "lon", "Year", "Month"))
df <- df %>% drop_na()

rm(species, env)


#### --- Prep the model --- ####

# --- Select species data

Y <- df %>%
  dplyr::select(6:162) # just species columns
Y = as.matrix(Y)

# --- Select environmental data

XData = data.frame(station = as.factor(df$Station), 
                   year = as.factor(df$Year), month = as.factor(df$Month),
                   o2 = df$o2, temp = df$temp,
                   ph = df$ph, depth = df$depth)


# --- Coordinate data

xy = as.matrix(cbind(df$lon, df$lat))
rownames(xy) = as.character(df$Station)
colnames(xy) = c("x-coordinate", "y-coordinate")
xy <- unique(xy)


#### --- Design model --- ####

# --- Study design

studyDesign = data.frame(station = XData$station, year = XData$year)


# --- Random effect structure (hierarchical study design)

rL.year = HmscRandomLevel(units = levels(studyDesign$year))
rL.station = HmscRandomLevel(units = levels(studyDesign$station))


# --- Regression model for environmental covariates

XFormula = ~ o2 + temp + ph + depth


#### --- Construct the model --- ####

if (arguments$model_type == "full"){
  model = Hmsc(Y = Y, XData = XData,
               XFormula = XFormula,
               distr = "probit", # because PA
               studyDesign = studyDesign,
               ranLevels = list(station = rL.station,
                                year = rL.year))
} else if (arguments$model_type == "environmental"){
  model = Hmsc(Y = Y, XData = XData, 
               XFormula = XFormula,
               distr = "probit")
} else if (arguments$model_type == "spatial"){
  model = Hmsc(Y = Y, XData = XData, 
               XFormula = ~1, # this is the difference
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(year = rL.year))
} else {stop(paste("unrecognised model type:", arguments$model_type))}


# Run model
sampler = sampleMcmc(model,
                     initPar = "fixed effects",
                     engine="HPC",
                     thin = arguments$n_thins, 
                     samples = arguments$n_samples, 
                     transient = arguments$n_transients,
                     nChains = arguments$n_chains
                     )

init_file_path = file.path(arguments$output_dir, arguments$sampler_rds)
saveRDS(to_json(sampler), file = init_file_path)

save(model, file = file.path(arguments$output_dir, arguments$model_rds))