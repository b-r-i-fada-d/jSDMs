library("argparse")

build.s01.argparser <- function(){
  parser <- ArgumentParser(description='HMSC-HPC Pipeline')
  parser$add_argument("-i", "--model-rds", type = "character", required=TRUE)
  parser$add_argument("-g", "--gibbs-samples-prefix", type = "character", required=TRUE)
  parser$add_argument("-o", "--output-dir", type = "character", required=TRUE)
  parser$add_argument("-m", "--model-type", type = "character", required=TRUE)
  
  parser$add_argument("-s", "--n-samples", type = "integer", required=TRUE)
  parser$add_argument("-c", "--n-chains", type = "integer", required=TRUE)
  parser$add_argument("-t", "--n-thins", type = "integer", required=TRUE)
  parser$add_argument("-r", "--n-transients", type = "integer", required=TRUE)
  
  parser$add_argument("-f", "--species-data", type = "character", required=FALSE)
  parser$add_argument("-e", "--env-data", type = "character", required=FALSE)
  parser$add_argument("-p", "--pa-data", type = "character", required=FALSE)
  
  # parser$add_argument("-x", "--seed", type = "integer", required=FALSE)
  parser$add_argument("-v", "--verbosity", type = "integer", required=FALSE)
  
  return(parser)
}

build.s02.argparser <- function(){
  parser <- ArgumentParser(description='HMSC-HPC Pipeline')
  parser$add_argument("-i", "--model-rds", type = "character", required=TRUE)
  # parser$add_argument("-g", "--sampler-rds", type = "character", required=TRUE)
  parser$add_argument("-o", "--output-dir", type = "character", required=TRUE)
  parser$add_argument("-m", "--model-type", type = "character", required=TRUE)
  
  parser$add_argument("-s", "--n-samples", type = "integer", required=TRUE)
  parser$add_argument("-c", "--n-chains", type = "integer", required=TRUE)
  parser$add_argument("-t", "--n-thins", type = "integer", required=TRUE)
  parser$add_argument("-r", "--n-transients", type = "integer", required=TRUE)
  
  return(parser)
}
