#!/bin/bash

#SBATCH --jobname=HMSC-HPC.Run
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --gpus-per-node=1
#SBATCH --time=00:20:00
#SBATCH --partition=gpu
#SBATCH --output=%x.%j.%a.out
#SBATCH --error=%x.%j.%a.err


set -euo pipefail

help () {
	cat <<- EOF

	${0} [-v <verbosity>] -g <sampler.rds> -o <output_dir> -s <n_samples> -t <n_thins> -r <n_transient>

	-g <sampler.rds>     Initialized Gibbs Sampler RDS file produced by R-HMSC. (path)
	-o <output_dir>      Output directory. Will be created if non-existant. (path)

	-s <n_samples>       Number of samples to generate. (int)
	-t <n_thins>         Number of thins. (int)
	-r <n_transient>     Number of transients. (int)
	
	-v <verbosity>       Verbosity level. (int)
	
	EOF
}

verbosity=0

while getopts ":g:o:s:t:r:v:h" arg; do
	case "${arg}" in
		g)
			sampler_rds="${OPTARG}"
			;;
		o)
			output_dir="${OPTARG}"
			;;
		s)
			n_samples="${OPTARG}"
			;;
		t)
			n_thins="${OPTARG}"
			;;
		r)
			n_transients="${OPTARG}"
			;;
		v)
			verbosity="${OPTARG}"
			;;
		h)
			help
			exit 0
			;;
		:)
			echo "Missing required argument for option: -${OPTARG}"
			help
			exit 1
			;;
		?)
			echo "Invalid option: -${OPTARG}"
			help
			exit 1
			;;
	esac
done

if [ ${OPTIND} -eq 1 ]; then help; exit 0; fi

task_id="${SLURM_ARRAY_TASK_ID:-1}"

python -m hmsc.run_gibbs_sampler \
	--input "${sampler_rds}" \
	--output "${output_dir}/${sampler_rds%.rds}_chain_${task_id}.json.gz" \
	--samples "${n_samples}" \
	--transient "${n_transients}" \
	--thin "${n_thins}" \
	--verbose "${verbosity}"
