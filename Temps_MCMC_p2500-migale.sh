#!/bin/bash
#$ -S /bin/bash
#Utilisateur  mnaveau
#$ -M marion.naveau@inrae.fr
#$ -m bea
#$ -cwd

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libopenblas.so"

R CMD BATCH Temps_MCMC_p2500.R Temps_MCMC_p2500.Rout
