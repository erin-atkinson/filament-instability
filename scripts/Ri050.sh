#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=filament-instability-Ri050
#SBATCH --output=../../scratch/logs/filament-instability.Ri050.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0

cd ~/Old/filament-instability

julia -t 32 -- src/simulation.jl default parameters/Ri050.txt 1 25 10 ../../scratch/filament-instability/Ri050