#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Prep
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --partition=largemem
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")
parser <- build.s02.argparser()
arguments <- parser$parse_args()

library(Hmsc)
library(tidyverse)
library(ggplot2)

# Read environmental data
grid <- read_csv(arguments$PA_data)

grid <- grid %>% select(-Year)

# Coordinates and environmental predictors
xy.grid <- as.matrix(cbind(grid$lon, grid$lat))
XData.grid <- grid %>% select(temp = BO2_templtmax_bdmean, 
                              o2 = BO_dissox, ph = BO_ph, 
                              depth = BO_bathymean, sal = BO_salinity)

# # Prepare gradient (if model has spatial effects)
# Gradient <- prepareGradient(model, XDataNew = XData.grid, sDataNew = list(station = xy.grid))

# Posterior predictive distribution
nParallel <- 2
predY <- predict(model, expected = TRUE, nParallel = nParallel)
EpredY <- Reduce("+", predY) / length(predY)

# Save predictions
save(EpredY, file = file.path("predictions.RData"))
# 
# # Post-process predictions
# DI <- EpredY[, 56]  # Example species
# S  <- rowSums(EpredY)  # Species richness
# CWM <- (EpredY %*% model$Tr) / matrix(rep(S, model$nt), ncol = model$nt)  # Community-weighted mean traits
# 
# # Prepare data frame for mapping
# mapData <- data.frame(
#   x = grid$lon,
#   y = grid$lat,
#   temp = XData.grid$temp,
#   o2 = XData.grid$o2,
#   ph = XData.grid$ph,
#   depth = XData.grid$depth,
#   DI = DI,
#   S = S,
#   CWM = CWM
# )
# 
# # Example clustering
# RCP10 <- kmeans(EpredY, 10)
# mapData$RCP10.cluster <- as.factor(RCP10$cluster)
# 
# # Create output folder for maps
# dir.create("maps", showWarnings = FALSE)
# 
# # Function to save a ggplot automatically
# save_plot <- function(plot_obj, filename) {
#   ggsave(filename = file.path("maps", filename),
#          plot = plot_obj, width = 6, height = 5, dpi = 300)
# }
# 
# # List of plots to save
# plots <- list(
#   ggplot(mapData, aes(x = x, y = y, color = temp)) +
#     geom_point(size = 1) + ggtitle("Temperature") + scale_color_gradient(low = "blue", high = "red") + coord_equal(),
#   
#   ggplot(mapData, aes(x = x, y = y, color = o2)) +
#     geom_point(size = 1) + ggtitle("Oxygen") + scale_color_gradient(low = "blue", high = "red") + coord_equal(),
#   
#   ggplot(mapData, aes(x = x, y = y, color = DI)) +
#     geom_point(size = 1) + ggtitle(expression(italic("Dipturus intermedius"))) + scale_color_gradient(low = "blue", high = "red") + coord_equal(),
#   
#   ggplot(mapData, aes(x = x, y = y, color = S)) +
#     geom_point(size = 1) + ggtitle("Species richness") + scale_color_gradient(low = "blue", high = "red") + coord_equal(),
#   
#   ggplot(mapData, aes(x = x, y = y, color = CWM)) +
#     geom_point(size = 1) + ggtitle("Community-weighted mean traits") + scale_color_gradient(low = "blue", high = "red") + coord_equal(),
#   
#   ggplot(mapData, aes(x = x, y = y, color = RCP10.cluster)) +
#     geom_point(size = 1) + ggtitle("RCP10 clusters") + scale_color_discrete() + coord_equal()
# )
# 
# # Filenames for saved plots
# filenames <- c("map_temp.png", "map_o2.png", "map_DI.png", "map_richness.png", "map_CWM.png", "map_RCP10.png")
# 
# # Save all plots
# for (i in seq_along(plots)) {
#   save_plot(plots[[i]], filenames[i])
# }
