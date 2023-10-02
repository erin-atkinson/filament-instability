#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=1:00:00
#SBATCH --job-name=filament-instability-default-test
#SBATCH --output=../scratch/logs/filament-instability/default_test.txt
module load cuda/11.0.3

cd ~/filament-instability

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src/simulation.jl default parameters/default.txt 1 10 100 ../scratch/filament-instability/default