#include<stdio.h>
#include<iostream>
#include <cstdlib>
#include<limits.h>
#include<algorithm>

using namespace std;

#define maxVertices   8192
#define INF           INT_MAX-1

float dist[maxVertices * maxVertices];
float *device_matrix;
float *result_matrix;
int vertices;
int tilesize[2];

void init(int n)
{
        for(int i=0;i<n;i++)
        {
                for(int j=0;j<n;j++)
                {
                        if(i==j)
                        {
                                dist[i * vertices + j] = 0;
                        }
                        else
                        {
                                dist[i * vertices + j] = INF;
                        }
                }
        }
}

void printarray(int vertices)
{
	return;

    for(int i = 0 ; i < vertices; i++ )
        {
                cout << "\n";
                for(int j = 0 ; j< vertices ;j++ )
                        cout << dist[i * vertices + j] << " " ;
        }

	printf("\n");
}

__global__
void FloydWarshall(int k, int idelta, int jdelta, float *matrix, int n, int na)
{
    int col = jdelta + blockIdx.x * blockDim.x + threadIdx.x; /* This thread’s matrix column */

    if(col >= n)
        return;

    int arrayIndex = n * (idelta + blockIdx.y) + col;

    __shared__ long trkc; /* this row, kth column */

    if(threadIdx.x == 0)
        trkc = matrix[na * (idelta + blockIdx.y) + k];

    __syncthreads();

    int tckr = matrix[k*na + col]; /* this column, kth row */

    int betterMaybe = trkc + tckr;

    if(betterMaybe < matrix[arrayIndex])
        matrix[arrayIndex] = betterMaybe;
}

void F_loop_FW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n)
{       
    int threads_per_block   = 1024;

    dim3 blocks_per_grid((n + threads_per_block - 1) /
                                threads_per_block, n);

        for(int via = Uj; via < Uj + n; via++)
        {
        FloydWarshall<<<blocks_per_grid, threads_per_block>>>(via, Xi, Xj, device_matrix, 
                     n, vertices);           
    }
}
/*
void F_loop_FW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n)
{       
	for(int via = Uj; via < Uj + n; via++)
	{
	for(int from = Xi; from < Xi + n; from++)
        {
                for(int to = Xj; to < Xj + n ; to++)
                {
                        if(from!=to && from!=via && to!=via)
			{
					dist[from * vertices + to] = min(dist[from * vertices + to],
				dist[from * vertices + via]+dist[via * vertices + to]);
			}
                        
                }
        }
   }

	printarray(vertices);
}
*/
void DFW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}


void BFW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    BFW(Xi + p, Xj + ip , Ui + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}

void CFW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    CFW(Xi + ip, Xj + p , Ui + ip, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}

void AFW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n, int d) {

    int r = tilesize[d];

	if (n < r)
		F_loop_FW(Xi, Xj, Ui, Uj, Vi, Vj, n);

	else {
        for (int k = 0; k < r; k++) {
            int p = k * (n/r);

            AFW(Xi + p, Xj + p, Ui + p, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    BFW(Xi + p, Xj + ip , Ui + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

            }

            for (int j = 0; j < r; j++) {
                int ip = j * (n/r);

                if (j != k)
                    CFW(Xi + ip, Xj + p , Ui + ip, Uj + p, Vi + p, Vj + p, n/r, d + 1);

            }

            for (int i = 0; i < r; i++)
               for (int j = 0; j < r; j++) {
                   int ip = i * (n/r);
                   int jp = j * (n/r);

                   if (i != k && j != k)
                       DFW(Xi + ip, Xj + jp, Ui + ip, Uj + p, Vi + p, Vj + jp, n/r, d + 1);
                }

        }
	}
}

int main(int argc, char *argv[])
{      
 	char *arg_vertices = getenv("N_VERTICES");
	vertices = atoi(arg_vertices);
	
    tilesize[0] = 2;
    tilesize[1] = INF;

    init(vertices);

	for(int i = 0 ; i < vertices ; i++ )
	{
		for(int j = 0 ; j< vertices; j++ )       
		{
			if( i == j )
				dist[i * vertices + j] = 0;
			else {
				int num = i + j;

				if (num % 3 == 0)
					 dist[i * vertices + j] = num / 2;
				else if (num % 2 == 0)
					 dist[i * vertices + j] = num * 2;
				else
					 dist[i * vertices + j] = num;
			}
		}
	}	

    size_t tot = vertices * vertices * sizeof(float);
    device_matrix = NULL;
    cudaMalloc((float **)&device_matrix, tot);

    cudaMemcpy(device_matrix, dist, tot, cudaMemcpyHostToDevice);


	AFW(0, 0, 0, 0, 0, 0, vertices, 0);

    result_matrix =(float *)malloc( vertices * vertices *
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
