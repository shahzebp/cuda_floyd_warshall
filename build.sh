module load cuda
nvcc -arch=compute_35 -code=sm_35 AloopFW.cu
rm Test*
