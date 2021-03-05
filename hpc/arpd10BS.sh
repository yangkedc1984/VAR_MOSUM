#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=72:00:00
module load languages/R-4.0.3-gcc9.1.0
export WORK_DIR=$HOME/sims
cd $WORK_DIR
Rscript arpd10BS.R

