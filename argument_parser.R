library("argparse")

parser <- ArgumentParser(description='HMSC-HPC Pipeline')
parser$add_argument("-i", "--model-rds", type = "character", required=TRUE)
parser$add_argument("-o", "--output-dir", type = "character", required=TRUE)
parser$add_argument("-m", "--model", type = "character", required=TRUE)

parser$add_argument("-s", "--n-samples", type = "integer", required=TRUE)
parser$add_argument("-c", "--n-chains", type = "integer", required=TRUE)
parser$add_argument("-t", "--n-thins", type = "integer", required=TRUE)
parser$add_argument("-r", "--n-transients", type = "integer", required=TRUE)

parser$add_argument("-f", "--species-data", type = "character", required=TRUE)
parser$add_argument("-e", "--env-data", type = "character", required=TRUE)

# parser$add_argument("-x", "--seed", type = "integer", required=FALSE)
parser$add_argument("-v", "--verbosity", type = "integer", required=FALSE)

arguments <- parser$parse_args()
print(arguments)