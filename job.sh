#!/bin/bash
#SBATCH -J Test        # Job Name
#SBATCH -o Test.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 16           # Total number of  tasks requested
#SBATCH -p gpu  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) - 1.5 hours

N_VERTICES=512 ./a.out

#N_VERTICES=$1 ./fwr1 > a_fwr1_output
#N_VERTICES=$1 ./fwr2 > a_fwr2_output
#N_VERTICES=$1 ./fwr3 > a_fwr3_output
