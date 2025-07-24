#!/usr/bin/env Rscript

#SBATCH --job-name=HMSC-HPC.Post
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source("argument_parser.R")

parser <- build.s02.argparser()
arguments <- parser$parse_args()

library("Hmsc")
library("tidyverse")
library("ggplot2")
library("vioplot")
library("colorspace")
library("jsonify")

samples_list = vector("list", arguments$n_chains) # empty list with slots
for(i in 0:(arguments$n_chains-1)){
  chain_file <- arguments$gibbs_samples_prefix
  chain_file <- file.path(paste0(chain_file, "_chain_", i, ".rds"))
  samples_list[[i+1]] = from_json(readRDS(file = chain_file)[[1]])[[1]]
}

model <- readRDS(file.path(arguments$model_rds))

model = importPosteriorFromHPC(model,
                                  samples_list, 
                                  arguments$n_samples, 
                                  arguments$n_thins, 
                                  arguments$n_transients)

fitted_filename <- gsub("\\.rds$", "_fitted.rds", arguments$model_rds)
saveRDS(model, file = file.path(fitted_filename))