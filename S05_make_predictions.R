#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Prep
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=20:00:00
#SBATCH --partition=largemem
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")
parser <- build.s03.argparser()
arguments <- parser$parse_args()

library(Hmsc)
library(tidyverse)
library(ggplot2)
library(parallel)

# Read environmental data
grid <- read_csv(arguments$env_data)
model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"

grid <- grid %>% drop_na()

nParallel <- 4

grid <- grid %>% dplyr::select(-geometry)
grid <- grid %>% rename(temp = SBT, month = Month, year = Year)

# xy.grid = as.matrix(cbind(grid$lon,
#                           grid$lat))
# 
# XData.grid <- data.frame(ph = grid$ph, 
#                          depth = grid$depth,
#                          o2 = grid$o2,
#                          temp = grid$temp,
#                          month = grid$month,
#                          year = grid$year,
#                          station = grid$station,
#                          stringsAsFactors = TRUE
# )
# 
# studyDesign.grid <- data.frame(
#   site = as.factor(1:nrow(XData.grid))
# )
# 
# 
# sDataNew <- list(station = XData.grid$station)
# 
# # Then prepare the gradient
# Gradient <- prepareGradient(
#   model,
#   XDataNew = XData.grid,
#   sDataNew = sDataNew
# )
# 
# #  --- Predict
# 
# # Use posterior mean prediction (ignoring MCMC uncertainty)
# predY <- predict(
#   model,
#   Gradient = Gradient,
#   expected = TRUE,
#   nParallel = 2
# )
# 
# # If predY is a list (from multiple chains/samples), average to get posterior mean
# EpredY <- Reduce("+", predY) / length(predY)
# 
# # Convert to long format for plotting or further analysis
# pred_df <- cbind(grid, EpredY) %>%
#   pivot_longer(
#     cols = starts_with("V"),  # adjust if species names are different
#     names_to = "species",
#     values_to = "predicted"
#   )
# 
# 
# write.csv(pred_df, "spatial_predictions.csv", row.names = FALSE)

# 
# new 16.09.2025
xy = as.matrix(cbind(grid$lon,
                     grid$lat))
xy <- unique(xy)

XData <- data.frame(ph = as.factor(grid$ph),
                    depth = as.factor(grid$depth),
                    o2 = as.factor(grid$o2),
                    temp = as.factor(grid$SBT),
                    month = as.factor(grid$Month),
                    year = as.factor(grid$Year),
                    # station = as.factor(grid$station),
                    stringsAsFactors = TRUE
)
# 
# studyDesign.grid <- data.frame(
#   station = as.factor(1:nrow(XData.grid))
# )


studyDesign <- data.frame(station = as.factor(grid$station),
                          year = as.factor(grid$year))


rL.station = HmscRandomLevel(sData = xy, sMethod = "NNGP")
rL.year = HmscRandomLevel(units = levels(studyDesign$year))

ranLevels = list(station = rL.station,
                 year = rL.year)


# old
# # Coordinates and environmental predictors
# xy.grid <- as.matrix(cbind(grid$lon, grid$lat))
# XData.grid <- grid %>% select(temp = temp, 
#                               o2 = o2, ph = ph, 
#                               depth = depth)
# # Prepare gradient (if model has spatial effects)
# Gradient <- prepareGradient(model, XDataNew = XData.grid, sDataNew = list(station = xy.grid))

# # Posterior predictive distribution
# 
# # old
# predY <- predict(model, expected = TRUE, nParallel = nParallel)
# EpredY <- Reduce("+", predY) / length(predY)
# 
# # Save predictions
# save(predY, file = file.path("predictions_raw.RData"))
# save(EpredY, file = file.path("predictions.RData"))

# new
#predY = predict(model, predictEtaMean = TRUE, expected = TRUE) # old
predY <- predict(model,
                 XData = XData,
                 studyDesign = studyDesign,
                 ranLevels = ranLevels,
                 predictEtaMean = TRUE,
                 expected = TRUE
)

EpredY.grid <- Reduce("+", predY.grid) / length(predY.grid)

save(EpredY.grid, file = "results/predictions_grid.RData")





#### --- Batching --- ####
# 
# batch_size <- 1000
# n_cores <- 4
# n_sites <- nrow(XData.grid)
# batches <- split(1:n_sites, ceiling(seq_along(1:n_sites)/batch_size))
# 
# predict_batch <- function(i) {
#   idx <- batches[[i]]
# 
#   cat("Running batch", i, "of", length(batches), "\n")
# 
#   predY.batch <- predict(
#     object = model,
#     XData = XData.grid[idx, , drop = FALSE],
#     studyDesign = studyDesign.grid[idx, , drop = FALSE],
#     ranLevels = list(),
#     expected = TRUE,
#     predictEtaMean = TRUE
#   )
# 
#   EpredY.batch <- Reduce("+", predY.batch) / length(predY.batch)
# 
#   batchDF <- data.frame(
#     site = idx,
#     lon = grid$lon[idx],
#     lat = grid$lat[idx],
#     EpredY.batch
#   )
# 
#   save(batchDF, file = paste0("results/predictions_batch_", i, ".RData"))
#   return(NULL)
# }
# 
# mclapply(seq_along(batches), predict_batch, mc.cores = n_cores)
# cat("All batches submitted.\n")
# 
# files <- list.files("results", pattern = "predictions_batch_.*\\.RData", full.names = TRUE)
# 
# batch_list <- lapply(files, function(f) {
#   load(f)  # loads batchDF
#   batchDF
# })
# 
# mapData <- do.call(rbind, batch_list)
# mapData <- mapData[order(mapData$site), ]
# mapData$site <- NULL
# 
# save(mapData, file = "results/predictions_grid_combined.RData")

