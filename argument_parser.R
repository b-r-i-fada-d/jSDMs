library(argparse)

parser <- ArgumentParser(description='HMSC')
parser$add_argument("-s", "--n_samples", type = "integer", required=TRUE)
parser$add_argument("-c", "--n_chains", type = "integer", required=TRUE)
parser$add_argument("-t", "--n_thins", type = "integer", required=TRUE)
parser$add_argument("-r", "--n_transients", type = "integer", required=TRUE)
parser$add_argument("-m", "--model", type = "character", required=TRUE)

arguments <- parser$parse_args()
