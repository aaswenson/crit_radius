#!/bin/sh
#SBATCH --partition=cnerg		# default "univ", if not specified
#SBATCH --nodes=1			# require 2 nodes
#SBATCH --ntasks-per-node=20            # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=4000		# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


export PATH=/home/group/dagmc/MCNP/bin:${PATH}
export DATAPATH=/home/group/dagmc/MCNP/MCNP_DATA

fuel='UO2'
cool='CO2'
matr='None'
frac='0.5'

mkdir ./${fuel}_${cool}_test${frac}
cp base_input.txt *py ./${fuel}_${cool}_test${frac}
cd ./${fuel}_${cool}_test${frac}

/home/aaswenson/python/bin/python3 critical_radius.py ${cool} ${fuel} Inconel-718 ${matr} ${frac}
