#include<stdio.h>
#include<iostream>
#include <cstdlib>
#include<limits.h>
#include<algorithm>

#include <cuda_runtime.h>

using namespace std;

#define maxVertices   8192
#define INF           INT_MAX-1

int vertices;
int tilesize[2];

void init(float *dist, int n)
{
        for(int i=0;i<n;i++)
        {
                for(int j=0;j<n;j++)
                {
                        if(i==j)
                        {
                                //dist[i][j] = 0;
                        }
                        else
                        {
                         //       dist[i][j] = INF;
                        }
                }
        }
}

__global__
void FloydWarshall(int k, int idelta, int jdelta, float *matrix, int n)
{
    int col = jdelta + blockIdx.x * blockDim.x + threadIdx.x; /* This thread’s matrix column */

    if(col >= n)
        return;

    int arrayIndex = n * (idelta + blockIdx.y) + col;

    __shared__ long trkc; /* this row, kth column */

    if(threadIdx.x == 0)
        trkc = matrix[n * (idelta + blockIdx.y) + k];

    __syncthreads();

    int tckr = matrix[k*n + col]; /* this column, kth row */

    int betterMaybe = trkc + tckr;

    if(betterMaybe < matrix[arrayIndex])
        matrix[arrayIndex] = betterMaybe;
}

/*
void F_loop_FW(float * device_matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n)
{       
    int threads_per_block   = 1024;

    dim3 blocks_per_grid((n + threads_per_block - 1) /
                                threads_per_block, vertices);

	for(int via = Uj; via < Uj + n; via++)
	{
        FloydWarshall<<<blocks_per_grid, threads_per_block>>>(via, Xi, Xj, device_matrix, 
                     vertices);           
    }
} */

void F_loop_FW(float *matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n)
{
    cout << "hitting " << endl;
    for(int via = Uj; via < Uj + n; via++)
    {
    for(int from = Xi; from < Xi + n; from++)
        {
                for(int to = Xj; to < Xj + n ; to++)
                {
                        if(from!=to && from!=via && to!=via)
            {
                 matrix[from * n + to] = min(matrix[from * n + to],
                            matrix[from * n + via] + matrix[via * n + to]);
            }
                 
                }
        }
   }
}

void DFW(float *matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(matrix, Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(matrix, Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}


void BFW(float *matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(matrix, Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    BFW(matrix, Xi + p, Xj + ip , Ui + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(matrix, Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}

void CFW(float *matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(matrix, Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    CFW(matrix, Xi + ip, Xj + p , Ui + ip, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(matrix, Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}

void AFW(float *matrix, int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(matrix, Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            AFW(matrix, Xi + p, Xj + p, Ui + p, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    BFW(matrix, Xi + p, Xj + ip , Ui + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

            }

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    CFW(matrix, Xi + ip, Xj + p , Ui + ip, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(matrix, Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
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

    tilesize[0] = 2;
    tilesize[1] = INF;
    //init(host_matrix, vertices);

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

    int threads_per_block   = 1024;
    dim3 blocks_per_grid((vertices + threads_per_block - 1) /
                                threads_per_block, vertices);

    AFW(host_matrix, 0, 0, 0, 0, 0, 0, vertices, 0);

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
