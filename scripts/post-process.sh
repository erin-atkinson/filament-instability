#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=7:00:00
#SBATCH --job-name=filament-instability-post-process
#SBATCH --output=../../scratch/logs/filament-instability.post-process.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara

# Path to ramdisk to avoid too many writes on parallel filesystem
mkdir /dev/shm/filament-instability-pp
export RAM=/dev/shm/filament-instability-pp

export PP_COMMAND=julia -t 40 -- src/analysis/post-process.jl

cd ~/Old/filament-instability

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../../scratch/filament-instability/Ri000

$PP_COMMAND $SIM_OUTPUT_FOLDER PSI_notnice $RAM

export SIM_OUTPUT_FOLDER=../../scratch/filament-instability/Ri010

$PP_COMMAND $SIM_OUTPUT_FOLDER PSI_nice $RAM
$PP_COMMAND $SIM_OUTPUT_FOLDER PSI_notnice $RAM

export SIM_OUTPUT_FOLDER=../../scratch/filament-instability/Ri020

$PP_COMMAND $SIM_OUTPUT_FOLDER PSI_nice $RAM
$PP_COMMAND $SIM_OUTPUT_FOLDER PSI_notnice $RAM
