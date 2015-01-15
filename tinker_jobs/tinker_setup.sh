#!/bin/bash

# This script sets up the Tinker directories for Mandor.
# It should be called by Settings.java.  If multiple jobs are being run,
# you should specify a node-local directory like /dev/shm.
#
# Usage: ./tinker_setup.sh /Users/ekwan/mandor2/tinker_jobs /Users/ekwan/mandor2/amino_acids /Users/ekwan/mandor2/tinker_jobs
#
# The first argument is the directory to setup the jobs into.
# The second argument is where to copy the forcefield files from.
# Typically this is amino_acids in the mandor working directory.
# The third argument is where to copy the tinker scripts from.
# Typically, this is tinker_jobs.
#
# The script will delete directories if they already exist in this directory
# and copy necessary files into them.

rm -rf $1/tinker_analysis_jobs
rm -rf $1/tinker_minimization_jobs

mkdir $1/tinker_analysis_jobs
mkdir $1/tinker_minimization_jobs

cp $2/amoebapro13.prm $1/tinker_analysis_jobs
cp $2/amoebapro13.prm $1/tinker_minimization_jobs

cp $2/oplsaal.prm $1/tinker_analysis_jobs
cp $2/oplsaal.prm $1/tinker_minimization_jobs

cp $3/run_tinker_analysis.sh $1/tinker_analysis_jobs
cp $3/run_tinker_minimization.sh $1/tinker_minimization_jobs
