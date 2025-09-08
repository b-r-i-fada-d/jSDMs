#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Prep
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --partition=largemem
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")
parser <- build.s03.argparser()
arguments <- parser$parse_args()

library(Hmsc)
library(tidyverse)
library(ggplot2)

# Read environmental data
grid <- read_csv(arguments$env_data)
model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"

grid <- grid %>% drop_na()

# Coordinates and environmental predictors
# old
# xy.grid <- as.matrix(cbind(grid$lon, grid$lat))
# XData.grid <- grid %>% select(temp = temp, 
#                               o2 = o2, ph = ph, 
#                               depth = depth)

# new
xy.grid = as.matrix(cbind(grid$lon,grid$lat))
XData.grid <- data.frame(ph = grid$ph, 
                         depth = grid$depth,
                         o2 = grid$o2,
                         temp = grid$temp,
                         month = grid$Month,
                         year = grid$Year,
                         stringsAsFactors = TRUE
)

studyDesign.grid <- data.frame(
  site = as.factor(1:nrow(XData.grid))
)

# # Prepare gradient (if model has spatial effects)
# Gradient <- prepareGradient(model, XDataNew = XData.grid, sDataNew = list(station = xy.grid))

# Posterior predictive distribution
nParallel <- 2

# old
# predY <- predict(model, expected = TRUE, nParallel = nParallel)
# EpredY <- Reduce("+", predY) / length(predY)
# 
# # Save predictions
# save(EpredY, file = file.path("predictions.RData"))

# new
# predY = predict(model, predictEtaMean = TRUE, expected = TRUE) # old
predY.grid <- predict(
  object = model,
  XData = XData.grid,
  studyDesign = studyDesign.grid,
  ranLevels = list(),
  predictEtaMean = TRUE,
  expected = TRUE
)

EpredY.grid <- Reduce("+", predY.grid) / length(predY.grid)


save(EpredY.grid, file = "results/predictions_grid.RData")
