#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=4:00:00
#SBATCH --job-name=filament-instability-post-process
#SBATCH --output=../../scratch/logs/filament-instability.post-process.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0

export SIM_OUTPUT_FOLDER=../../scratch/filament-instability/Ri010

cd ~/Old/filament-instability

julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER DFM
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER PV
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER v_balance