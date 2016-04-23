module load cuda
nvcc -arch=compute_35 -code=sm_35 $1 -o $2
rm Test*
