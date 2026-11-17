#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 8
#SBATCH --mem=50GB
#SBATCH -t 24:00:00
#SBATCH -J psrobot

psrobot="apptainer exec -B ${BIND_DIR}:${BIND_DIR} ${CONTAINER_DIR}/psrobot.sif"

$psrobot psRobot_tar -s $1 -t $2 -o $3
