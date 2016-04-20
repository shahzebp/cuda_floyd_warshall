#include<stdio.h>
#include<iostream>
#include <cstdlib>
#include <cilkview.h>
#include<limits.h>
#include<algorithm>
#include<cilk/cilk.h>

using namespace std;

#define maxVertices   8192
#define INF           INT_MAX-1

int dist[maxVertices][maxVertices];
int vertices;
int tilesize[2];

void init(int n)
{
        cilk_for(int i=0;i<n;i++)
        {
                cilk_for(int j=0;j<n;j++)
                {
                        if(i==j)
                        {
                                dist[i][j] = 0;
                        }
                        else
                        {
                                dist[i][j] = INF;
                        }
                }
        }
}

void F_loop_FW(int Xi, int Xj, int Ui, int Uj, int Vi, int Vj, int n)
{       
	for(int via = Uj; via < Uj + n; via++)
	{
	cilk_for(int from = Xi; from < Xi + n; from++)
        {
                cilk_for(int to = Xj; to < Xj + n ; to++)
                {
                        if(from!=to && from!=via && to!=via)
			{
				dist[from][to] = min(dist[from][to],dist[from][via]+dist[via][to]);
			}
                        
                }
        }
   }
}

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
                    BFW(Xi + p, Xj + ip , U + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

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
                    BFW(Xi + p, Xj + ip , U + p, Uj + p, Vi + p, Vj + ip, n/r, d + 1);

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
				dist[i][j] = 0;
			else {
				int num = i + j;

				if (num % 3 == 0)
					 dist[i][j] = num / 2;
				else if (num % 2 == 0)
					 dist[i][j] = num * 2;
				else
					 dist[i][j] = num;
			}
		}
	}	

	AFW(0, 0, 0, 0, 0, 0, vertices, 0);
	
    for(int i = 0 ; i < vertices; i++ ) 
	{
		cout << "\n";
		for(int j = 0 ; j< vertices ;j++ )
			cout << dist[i][j] << " " ;
	}

	return 0;
}
