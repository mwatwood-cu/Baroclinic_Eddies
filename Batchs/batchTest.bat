#!/bin/bash

#SBATCH --account=blanca-igg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --output=mesoscale-%j.out

echo "Starting"
../software/julia-1.3.0/bin/julia RunScript.jl
echo "Ending"
