#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=2:00:00
#SBATCH --job-name=filament-instability-cooling1
#SBATCH --output=../scratch/logs/filament-instability/cooling1.txt
module load cuda/11.0.3

cd ~/filament-instability

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src/simulation.jl secondarycirculation parameters/Turbulence/1.txt 30 0 100 ../scratch/filament-instability/cooling1