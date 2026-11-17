#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 2
#SBATCH --mem=20GB
#SBATCH -t 168:00:00

tapir="apptainer exec -B ${BIND_DIR}:${BIND_DIR} ${CONTAINER_DIR}/tapir.sif"

# it looks like this can't run in parallel. I try 4 cores and 8 cores, they run at the same speed (which is very slow)

# Usage hybrid_parser [options]
#    --mimic <1> search for target mimics (default 0)
#    --score <value> score cutoff (default 4)
#    --mfe <value> mfe ratio cutoff (default 0.7)

# we want to find target do I will not turn on mimic

$tapir tapir_hybrid $1 $2 | $tapir hybrid_parser > $3
