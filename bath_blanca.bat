#!/bin/bash

#SBATCH --account=blanca-igg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00
#SBATCH --output=mesoscale-%j.out

echo "Starting"
time ../software/julia-1.3.0/bin/julia RunScript_blanca.jl
echo "Ending"
