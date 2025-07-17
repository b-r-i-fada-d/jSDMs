#!/bin/bash

set -euo pipefail

help () {
	cat <<- EOF

	${0} [-v <verbosity>] -i <model.rds> -g <sampler.rds> -o <output_dir>
	-m <model_type> -f <species.csv> -e <env.csv>
	-c <n_chains> -s <n_samples> -t <n_thins> -r <n_transient>

	-i <model.rds>       Model RDS file to produce from R-HMSC. (path)
	-g <sampler.rds>     Initialized Gibbs Sampler RDS file to produce from R-HMSC. (path)
	-o <output_dir>      Output directory. Will be created if non-existant. (path)

	-c <n_chains>        Number of chains to run in parallel. (int)
	-s <n_samples>       Number of samples to generate. (int)
	-t <n_thins>         Number of thins. (int)
	-r <n_transient>     Number of transients. (int)

	-m <model_type>      Model type to run ("full", "spatial", or "environmental")
	-f <species.csv>     Species data file. (path)
	-e <env.csv>         Environmental data file. (path)

	-v <verbosity>       Verbosity level. (int)
	
	EOF
}

verbosity=0

while getopts ":i:g:o:m:c:s:t:r:f:e:v:h" arg; do
	case "${arg}" in
		i)
			model_rds="${OPTARG}"
			;;
		g)
			sampler_rds="${OPTARG}"
			;;
		o)
			output_dir="${OPTARG}"
			;;
		m)
			model_type="${OPTARG}"
			;;
		c)
			n_chains="${OPTARG}"
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
		f)
			species_data="${OPTARG}"
			;;
		e)
			env_data="${OPTARG}"
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


# 1) Submit the R script that prepares the model.
prep_job=$(sbatch --parsable \
	S01_fit_models.R \
		--model-rds "${model_rds}" \
		--sampler-rds "${sampler_rds}" \
		--output-dir "${output_dir}" \
		--n-samples "${n_samples}" \
		--n-chains "${n_chains}" \
		--n-thins "${n_thins}" \
		--n-transients "${n_transients}" \
		--model "${model_type}" \
		--species-data "${species_data}" \
		--env-data "${env_data}" \
		--verbosity "${verbosity}"
	)

# 2) Submit a job for each python sampling chain. These will wait to run until the prep script finishes.
run_job=$(sbatch --parsable --dependency=afterok:"${prep_job}" --array=1-"${n_chains}" \
	hmsc-hpc.run.sh \
		-g "${output_dir}/${sampler_rds}" \
		-o "${output_dir}" \
		-s "${n_samples}" \
		-t "${n_thins}" \
		-r "${n_transients}" \
		-v "${verbosity}" \
	)

# 3) Submit the R script to analyze the results. This will wait to run until the chains are all finished.
sbatch --dependency=afterok:"${run_job}" hmsc-hpc.finalize.sh
