#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00
#SBATCH --partition=shas-testing

source /curc/sw/anaconda3/2019.03/bin/activate
conda activate julia
julia RunScript.jl