#!/bin/bash
#SBATCH -J Test        # Job Name
#SBATCH -o Test.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 16           # Total number of  tasks requested
#SBATCH -p development  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) - 1.5 hours

m=`echo "2^13" | bc`
for (( c = 16; c <= $m; c = c * 2))
do
        export N_VERTICES=$c
        echo -n "N_VERTICES = " $c " " >> $1
        ./Cloop_CPU >> $1
done
