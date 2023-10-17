#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=2:00:00
#SBATCH --job-name=filament-instability-Ro2Ri1Ek0
#SBATCH --output=../scratch/logs/filament-instability/Ro2Ri1Ek0.txt
module load cuda/11.0.3

cd ~/filament-instability

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src/simulation.jl secondarycirculation parameters/ParameterSweep/Ro2Ri1Ek0.txt 1 10 100 ../scratch/filament-instability/Ro2Ri1Ek0
