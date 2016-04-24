#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <algorithm>
#include <sys/time.h>
#include <cuda_runtime.h>

using namespace std;

#define INF           INT_MAX-1

__global__
void FloydWarshall(int k, int i, float *matrix, int n)
{
    int col = blockIdx.x * blockDim.x + threadIdx.x; /* This thread’s matrix column */

    if(col >= n)
        return;

    int arrayIndex = n * i + col; // Calculating D[i][j]

    float matrix_value_i_k = matrix[n * i + k];   // Calculating D[i][k]
    
    float matrix_value_k_j = matrix[k*n + col]; /* this column, kth row */
    
    matrix[arrayIndex] = fmin(matrix[arrayIndex], matrix_value_i_k + matrix_value_k_j);
}



int main(int argc, char *argv[])
{      
    char *arg_vertices = getenv("N_VERTICES");
    char *arg_threads_per_block = getenv("N_THREADS");
	
    size_t vertices = atoi(arg_vertices);
    int threads_per_block   = atoi(arg_threads_per_block);
   
    float *host_matrix =(float *)malloc( vertices * vertices *
                sizeof(float));
    
	for(int i = 0 ; i < vertices ; i++ ) {
		for(int j = 0 ; j< vertices; j++ ) {
            if( i == j )
                host_matrix[i * vertices + j] = 0;

            else {
				int num = i + j;

				if (num % 3 == 0)
					 host_matrix[i * vertices + j] = num / 2;
				else if (num % 2 == 0)
					 host_matrix[i * vertices + j] = num * 2;
				else
					 host_matrix[i * vertices + j] = num;
			}
		}
	}	

    
    size_t tot = vertices * vertices * sizeof(float);
    float *device_matrix = NULL;
    cudaMalloc((float **)&device_matrix, tot);

    cudaMemcpy(device_matrix, host_matrix, tot, cudaMemcpyHostToDevice);

    int blocks_per_grid = vertices + (threads_per_block - 1) /threads_per_block;
    struct timeval tvalBefore, tvalAfter;
    gettimeofday (&tvalBefore, NULL);
    for(int via = 0; via < vertices; via++) {
	for(int j = 0; j < vertices; j++){

	    FloydWarshall<<<blocks_per_grid, threads_per_block>>>(via, j, device_matrix, vertices);
        cudaThreadSynchronize();
	}
    }
    gettimeofday (&tvalAfter, NULL);
    printf("Time: %ld microseconds\n",
        ((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
        +tvalAfter.tv_usec) - tvalBefore.tv_usec
        );
    float *result_matrix =(float *)malloc( vertices * vertices *
                sizeof(float));
 
    cudaMemcpy(result_matrix, device_matrix, tot, cudaMemcpyDeviceToHost);
    return 0;
}
