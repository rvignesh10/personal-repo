#!/bin/sh
#SBATCH --ntasks=16
#SBATCH --nodes=4
#SBATCH --time=02:00:00
#SBATCH --output=op.txt

mpirun ./heat2d.bin -nx=1000 -tFinal=0.01 -debug=0 -saveMatlab=0 -option=2
