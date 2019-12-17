#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=shas
#SBATCH --output=mesoscale-%j.out

echo "Starting"
../software/julia-1.3.0/bin/julia RunScript.jl
echo "Ending"
