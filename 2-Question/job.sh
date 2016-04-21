#!/bin/bash
#SBATCH -J Test        # Job Name
#SBATCH -o Test.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 16           # Total number of  tasks requested
#SBATCH -p gpu  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) - 1.5 hours

N_THREADS=1024 N_VERTICES=256 ./a.out
