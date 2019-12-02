// ConsoleApplication37.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#define m_ptintf if (myrank==0)printf
#define L 1000
#define ITMAX 100

int i, j, k, it;
int ll, shift;
double(*A)[L];
double(*B)[L];
int main()
{
	MPI_Request req[4];
	int myrank, ranksize;
	int startrow, lastrow, nrow;
	MI_Status status[4];
	double t1, t2, time;
	MI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
	MPI_Barrier(MPI_COMM_WORLD);

    startrow = (myrank * L) / ranksize;
	lastrow = (((myrank + 1)*L) / ranksize) - 1;
	nrow = lastrow - startrow + 1;
	m_prntf("1 start\n");

	A = malloc((nrow + 2) * L * sizeof(double));
	B = malloc((nrow) * L * sizeof(double));
	for ( int i = 1; i <= nrow; i++)
        for( j = 0; j <= L-1; j++)
    {
        A[i][j] = 0;
        B[i-1][j] = 1. + startrow + i - 1 +j;
    }

    t1=MPI_Wtime();
    for (it = 1; it <= ITMAX; it++)
    {
        for (i = 1; i <= nrow; i++)
        {
            if(((i==1)&&(myrank==0))||(i==nrow)&&(myrank==ranksize-1))
                continue;
            for ( j = 1; j <= L - 2; j++ )
            {
               A[i][j] = B[i-1][j];

            }
        }

        if (myrank!=0)
            MPI_Irecv(&A[0][0],L, MPI_DOUBLE, myrank-1, 1235, MPI_COMM_WORLD, &req[0]);

        if (myrank!=ranksize-1)
            MPI_Isend(&A[nrow][0],L, MPI_DOUBLE, myrank + 1, 1235, MPI_COMM_WORLD, &req[2]);

        if (myrank!=ranksize-1)
            MPI_Irecv(&A[nrow+1][0],L, MPI_DOUBLE, myrank+1, 1236, MPI_COMM_WORLD, &req[3]);

        if (myrank!=0)
            MPI_Isend(&A[1][0],L, MPI_DOUBLE, myrank-1, 1236, MPI_COMM_WORLD, &req[1]);

            ll = 4; shift = 0;
            if (myrank==0){ll = 2; shift = 2;}
            if (myrank==ranksize-1){ll = 2;}
            MPI_Waitall(ll, &req[shift], &status[0]);

            for (i=1; i <= nrow; i++)
            {

                if(((i==1)&&(myrank=0))||((i=nrow)&&(myrank==ranksize-1)))
                    continue;

                for(j = 1; j <= L - 2; j++)
                    B[i-1][j] = (A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1])/4.;

            }
    }
    printf("%d: Time of task = %if\n", myrank, MPI_Wtime()-t1);
    MPI_Finalize();
    return 0;
}

