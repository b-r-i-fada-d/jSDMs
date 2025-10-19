#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Predict
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=4:00:00
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

nParallel <- 4

# Read environmental data
grid <- read_csv(arguments$env_data)
model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"

grid <- grid %>% drop_na()
grid <- grid %>% rename(temp = SBT, month = Month, year = Year)

grid <- grid %>% filter(year >= 2021)

grid <- grid %>% sample_n(10) # testing small subset

XData <- data.frame(station = as.factor(grid$station), #19.10.25 added in
                    lat = as.factor(grid$lat),
                    lon = as.factor(grid$lon),
                    year = as.factor(grid$year),
                    month = as.factor(grid$month),
                    o2 = grid$o2,
                    temp = grid$temp,
                    ph = grid$ph,
                    depth = grid$depth)

xy = as.matrix(cbind(XData$lon, XData$lat))
rownames(xy) = as.character(XData$station)
colnames(xy) = c("x-coordinate", "y-coordinate")
# xy <- unique(xy) 


##########################################################################
#### --- 13.10.2025 - trying all years --- ####
# # grid <- grid %>% filter(year == 2021)
# 
# grid1 <- grid
# grid1 <- grid1 %>% filter(month == 11)
# 
# XData <- data.frame(year = as.factor(grid1$year),
#                     month = as.factor(grid1$month),
#                     o2 = grid1$o2,
#                     temp = grid1$temp,
#                     ph = grid1$ph,
#                     depth = grid1$depth)
# # rownames(XData) <- c(1:nrow(XData))
# 
# xy = as.matrix(cbind(grid1$lon, grid1$lat))
# rownames(xy) = as.character(grid1$station)
# colnames(xy) = c("x-coordinate", "y-coordinate")
# xy <- unique(xy)


##########################################################################
# #### --- 13.10.2025 - try with gradient --- ####
# # --- 19.10.25 removing gradient
# # - gradient work sbut prediction takes too long
# 
# Gradient = prepareGradient(model, 
#                            XData = XData, 
#                            sData = list(station = xy))
# 
# start_time <- Sys.time()
# cat("Start time:", start_time, "\n")
# 
# predY <- Hmsc:::predict.Hmsc(model,
#                              XData = XData,
#                              Gradient = Gradient,
#                              expected = T,
#                              predictEtaMean = T,
#                              thin = 10)
# 
# end_time <- Sys.time()
# cat("End time:", end_time, "\n")
# cat("Elapsed time:", end_time - start_time, "\n")


################################################################
# #13.10.2025 commenting out to try with gradient
# #19.10.25 adding back in, removing gradientm removing year
# #13.10.2025 adding spatial & year rLs
# studyDesign <- data.frame(station = as.factor(XData$station),
#                           #station = xy,
#                           # year = as.factor(grid$year),
#                           stringsAsFactors = TRUE)
# 
# 
# rL.station = HmscRandomLevel(sData = xy, sMethod = "NNGP")
# # ranLevels = list(station = rL.station)
# #rL.year = HmscRandomLevel(units = levels(studyDesign$year))
# 
# 
# rownames(rL.station$s) <- XData$station #19.10.25
# 
# 
# ranLevels = list(station = rL.station#,
#                  #year = rL.year
# )

################################################################


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
# # new 16.09.2025
# xy = as.matrix(cbind(grid$lon,
#                      grid$lat))
# xy <- unique(xy)
# 
# XData <- data.frame(ph = as.factor(grid$ph),
#                     depth = as.factor(grid$depth),
#                     o2 = as.factor(grid$o2),
#                     temp = as.factor(grid$SBT),
#                     month = as.factor(grid$Month),
#                     year = as.factor(grid$Year),
#                     # station = as.factor(grid$station),
#                     stringsAsFactors = TRUE
# )
# 
# studyDesign.grid <- data.frame(
#   station = as.factor(1:nrow(XData.grid))
# )


# studyDesign <- data.frame(station = as.factor(grid$station),
#                           year = as.factor(grid$year))
# 
# 
# rL.station = HmscRandomLevel(sData = xy, sMethod = "NNGP")
# rL.year = HmscRandomLevel(units = levels(studyDesign$year))
# 
# ranLevels = list(station = rL.station,
#                  year = rL.year)


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

#############################################################################
# # TRY THIS 24.09.2025 - it worked but there was a typo while saving
#predY = predict(model, predictEtaMean = TRUE, expected = TRUE) # old


#############################################################################

# 13.10.2025 commenting out to try with grandient
# 19.10.25 trying without gradient
predY <- predict(model,
                 XData = XData,
                 studyDesign = studyDesign,
                 ranLevels = ranLevels#,
                 # predictEtaMean = TRUE, # 13.10 commented out
                 # expected = TRUE # 13.10 commented out
)

#############################################################################

EpredY <- Reduce("+", predY) / length(predY)

save(EpredY, file = file.path("results_station_randomlevel/predictions_station_randomlevel.RData"))
# #############################################################################

# 
# 
# 
# library(future.apply)
# library(abind)
# 
# # -----------------------------
# # Batching setup
# # -----------------------------
# batch_size <- 1000
# n_sites <- nrow(XData)
# batches <- split(seq_len(n_sites), ceiling(seq_len(n_sites) / batch_size))
# 
# cat("Total sites:", n_sites, "\n")
# cat("Batch size:", batch_size, "\n")
# cat("Number of batches:", length(batches), "\n\n")
# 
# # -----------------------------
# # Set up parallel backend
# # -----------------------------
# n_cores <- 4   # adjust to your machine
# plan(multisession, workers = n_cores)
# 
# # -----------------------------
# # Function to run prediction on one batch
# # -----------------------------
# predict_batch <- function(i) {
#   idx <- batches[[i]]
#   cat("Starting batch", i, "(", length(idx), "rows ) ...\n")
#   
#   predY.batch <- predict(
#     object = model,
#     XData = XData[idx, , drop = FALSE],
#     studyDesign = studyDesign[idx, , drop = FALSE],
#     ranLevels = ranLevels,
#     predictEtaMean = TRUE,
#     expected = TRUE
#   )
#   
#   # average across posterior samples (faster than Reduce("+", ...))
#   EpredY.batch <- apply(simplify2array(predY.batch), c(1, 2), mean)
#   
#   cat("Finished batch", i, "\n")
#   return(EpredY.batch)
# }
# 
# # -----------------------------
# # Run predictions in parallel
# # -----------------------------
# batch_results <- future_lapply(seq_along(batches), predict_batch,
#                                future.seed = TRUE)
# 
# cat("\nAll batches completed.\n")
# 
# # -----------------------------
# # Recombine batches
# # -----------------------------
# EpredY <- do.call(rbind, batch_results)
# 
# # -----------------------------
# # Save final result
# # -----------------------------
# save(EpredY, file = "results/predictions_noyr.RData")
# cat("Saved predictions to results/predictions_noyr.RData\n")


# # #### --- Batching --- ####
# 
# batch_size <- 1000
# n_cores <- 4
# n_sites <- nrow(XData)
# batches <- split(1:n_sites, ceiling(seq_along(1:n_sites)/batch_size))
# 
# predict_batch <- function(i) {
#   idx <- batches[[i]]
# 
#   cat("Running batch", i, "of", length(batches), "\n")
# 
#   # predY.batch <- predict(
#   #   object = model,
#   #   XData = XData[idx, , drop = FALSE],
#   #   studyDesign = studyDesign[idx, , drop = FALSE],
#   #   ranLevels = ranLevels,
#   #   expected = TRUE,
#   #   predictEtaMean = TRUE
#   # )
#   
#   predY.batch <- Hmsc::predict.Hmsc(model,  # removed one :
#                                XData = XData[idx, , drop = FALSE],
#                                Gradient = Gradient,
#                                expected = T, 
#                                predictEtaMean = T)
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
#   save(batchDF, file = paste0("results_spatiotemporal_randomlevels/predictions_batch_", i, ".RData"))
#   return(NULL)
# }
# 
# mclapply(seq_along(batches), predict_batch, mc.cores = n_cores)
# cat("All batches submitted.\n")
# 
# files <- list.files("results_spatiotemporal_randomlevels", pattern = "predictions_batch_.*\\.RData", full.names = TRUE)
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
# save(mapData, file = "results_spatiotemporal_randomlevels/predictions_grid_combined.RData")
# 
