#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <algorithm>

#include <cuda_runtime.h>

using namespace std;

#define INF           INT_MAX-1

void init(float *matrix, int n)
{
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            if(i==j)
            {
                matrix[i * n + j] = 0;
            }
            else
            {
                matrix[i * n + j] = INF;
            }
        }
    }
}

__global__
void FloydWarshall(int k, int i, float *matrix, int n)
{
    int col = blockIdx.x * blockDim.x + threadIdx.x; /* This thread’s matrix column */

    if(col >= n)
        return;

    int arrayIndex = n * i + col; // Calculating D[i][j]

    float trkc = matrix[n * i + k];   // Calculating D[i][k]
    
    float tckr = matrix[k*n + col]; /* this column, kth row */
    
    matrix[arrayIndex] = fmin(matrix[arrayIndex], trkc + tckr);
}



int main(int argc, char *argv[])
{      
    char *arg_vertices = getenv("N_VERTICES");
    char *arg_threads_per_block = getenv("N_THREADS");
	
    size_t vertices = atoi(arg_vertices);
    int threads_per_block   = atoi(arg_threads_per_block);
   
    float *host_matrix =(float *)malloc( vertices * vertices *
                sizeof(float));
    
	init(host_matrix, vertices);
    
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
