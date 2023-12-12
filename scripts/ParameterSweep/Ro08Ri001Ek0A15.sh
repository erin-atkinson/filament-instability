#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=filament-instability-Ro08Ri001Ek0A15
#SBATCH --output=../scratch/logs/filament-instability/Ro08Ri001Ek0A15.txt
module load cuda/11.0.3

cd ~/filament-instability

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src/simulation.jl secondarycirculation parameters/ParameterSweep/Ro08Ri001Ek0A15.txt 1 100 100 ../scratch/filament-instability/Ro08Ri001Ek0A15
