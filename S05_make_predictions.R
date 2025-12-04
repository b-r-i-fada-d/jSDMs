#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Predict
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=30:00:00
#SBATCH --partition=largemem
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R") 
parser <- build.s03.argparser()
arguments <- parser$parse_args()

library("Hmsc")
library("tidyverse")
library("ggplot2")
library("parallel")

nParallel <- 40

# Read environmental data
grid <- read_csv(arguments$env_data)
model <- readRDS(file.path(gsub("\\.rds$", "_fitted.rds", arguments$model_rds)))
output_prefix <- "model"

grid <- grid %>% drop_na()

xy = as.matrix(cbind(grid$lon, grid$lat))
rownames(xy) = as.character(grid$station)
colnames(xy) = c("x-coordinate", "y-coordinate")

XData.grid <- data.frame(year = as.factor(grid$year),
                    # month = as.factor(grid$month),
                    o2 = grid$o2,
                    temp = grid$temp,
                    ph = grid$ph,
                    depth = grid$depth,
                    stringsAsFactors = TRUE)
# Prepare gradient
Gradient = prepareGradient(model, 
                           XDataNew = XData.grid, 
                           sDataNew = list(station = xy))

# Compute the posterior predictive distribution
start_time <- Sys.time()
cat("Start time:", start_time, "\n")

predY = predict(model, 
                Gradient = Gradient, 
                predictEtaMean = TRUE,
                expected = TRUE, 
                nParallel=nParallel)

end_time <- Sys.time()
cat("End time:", end_time, "\n")
cat("Elapsed time:", end_time - start_time, "\n")

# predY is a list (from multiple chains/samples), average to get posterior mean
EpredY <- Reduce("+", predY) / length(predY)

save(predY, file = file.path("results_2014-2024_no_month_1000_bottom_env/Predictions_raw.RData"))
save(EpredY, file = file.path("results_2014-2024_no_month_1000_bottom_env/predictions.RData"))
