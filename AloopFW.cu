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
void FloydWarshall(float *matrix, int n)
{
        
	for(int via=0; via < n; via++) {
	    for(int from=0;from<n;from++) {
            for(int to=0;to<n;to++) {
                if(from!=to && from!=via && to!=via) {

				    matrix[from * n + to] = min(matrix[from * n + to], 
                            matrix[from * n + via] + matrix[via * n + to]);

			    }
                        
            }
        }
   	}
}


int main(int argc, char *argv[])
{      
    char *arg_vertices = getenv("N_VERTICES");
	
    size_t vertices = atoi(arg_vertices);
   
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

	//FloydWarshall(host_matrix, vertices);
	FloydWarshall<<<1, 1>>>(device_matrix, vertices);

    float *result_matrix =(float *)malloc( vertices * vertices *
                sizeof(float));
 
    cudaMemcpy(result_matrix, device_matrix, tot, cudaMemcpyDeviceToHost);
    
    for(int i = 0 ; i < vertices; i++ ) 
	{
		cout << "\n";
		for(int j = 0 ; j< vertices ;j++ )
			cout << result_matrix[i * vertices + j] << " " ;
	} 

	return 0;
}
