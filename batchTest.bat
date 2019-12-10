#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00
#SBATCH --partition=shas-testing
#SBATCH --output=mesoscale-out.out

module load gcc/8.2.0
module load python/3.6.5 cmake/3.14.1 
julia RunScript.jl
