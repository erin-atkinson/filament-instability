#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=1:00:00
#SBATCH --job-name=filament-instability-VaryingRi06
#SBATCH --output=../scratch/logs/filament-instability/VaryingRi06.txt
module load cuda/11.0.3

cd ~/filament-instability

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src/simulation.jl default parameters/VaryingRi/06.txt 1 10 100 ../scratch/filament-instability/VaryingRi06