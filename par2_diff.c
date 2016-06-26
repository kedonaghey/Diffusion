#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int rank, size;
int P, Q;
int lrows, lcols;

void InitGrid(double** mat)
{
  int i;
  if(rank%P == 0)
  {
    for(i= 0;i<lrows;i++)
	mat[i][1] = 20;
  }
  if(rank%P == P-1)
  {
    for(i=0;i<lrows ;i++)
	mat[i][lcols-2] = 30;
  }
  if(rank/P == 0)
  {
    for(i=0 ;i<lcols;i++)
	mat[1][i] = 50;
  }
  if(rank/P == Q-1)
  {
    for(i=0;i<lcols;i++)
	mat[lrows - 2][i] = 40;
  }
}

void PrintGrid(double** mat, int nrows, int ncols)
{
/*
int i, j, k, l; 
//  for(l=0;l<Q;l++) {
    for (i = 1; i < lrows -1; i++) {
//	for(k=0;k<P;k++) {
        for (j = 1; j < lcols -1; j++) {
            printf("%f ", mat[i][j]);
        }
  //      }
        printf("\n");
      }
//}
  */  

}
void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{
  int i, j, m;
  double **tmp;
  tmp = calloc(lrows, sizeof(double *));
  for (i = 0; i < lrows; i++)
    tmp[i] = calloc(lcols, sizeof(double));
  double wx = 0, wy = 0;
  //diffusitivity coefficient - can equal .002
  double k = .02;
  double d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0, d5 = 0.0;

  //weight in x direction
  wx = k * dt/(dx*dx);
  //weight in y direction
  wy = k * dt/(dy*dy);
  
  //stability
  // wx + wy must be less than .5
  if ( wx + wy > .5)
  {
    printf("Does not reach requirements for stability\n");
    exit(1);
  }

  //MPI
  MPI_Datatype MPI_CLM;
  MPI_Type_vector(lrows, 1, lcols, MPI_DOUBLE, &MPI_CLM);
  MPI_Type_commit(&MPI_CLM);
  MPI_Status stat[8];
  MPI_Request req[8];
  int u, d, l, r;
  int startrow, startcol, endrow, endcol;

  //sets when receiver/sender process doesnt exist
  u = rank - P;
  if ( u < 0) 
    u = MPI_PROC_NULL;
  d = rank + P;
  if (d >= size)
    d = MPI_PROC_NULL;
  l = rank - 1;
  if (rank%P == 0)
    l = MPI_PROC_NULL;
  r = rank + 1;
  if (rank%P == P - 1)
    r = MPI_PROC_NULL;



  if (rank%P == 0)
    startcol = 2;
  else
    startcol = 1;
  if (rank%P == P - 1)
    endcol = lcols - 2;
  else
    endcol = lcols - 1;
  if (rank/P == 0)
    startrow = 2;
  else
    startrow = 1;
  if (rank/P == Q - 1)
    endrow = lrows - 2;
  else
    endrow = lrows - 1;

  //solving edges
  if (rank/P != 0)
  {
    //not in top blocks
    i = 1;
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d1 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d1 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  if (rank/P != Q - 1)
  {
    //not in bottom blocks
    i = lrows - 2;
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d2 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d2 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  if (rank%P != 0)
  {
    //not in left blocks
    j = 1;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d3 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d3 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  if (rank%P != P -1)
  {
    //not in right blocks
    j = lcols - 2;
    for(i = startrow; i < endrow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d4 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d4 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }
  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[1], lcols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[0]);
  MPI_Irecv(mat2[lrows-1], lcols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[1]);
  MPI_Isend(mat2[lrows-2], lcols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[2]);
  MPI_Isend(mat2[2], lcols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[0][1], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[4]);
  MPI_Irecv(&mat2[0][lcols-1], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[5]);
  MPI_Isend(&mat2[0][2], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[6]);
  MPI_Isend(&mat2[0][lcols-2], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[7]);

  for(i = startrow; i < endrow; i++)
  {
    for(j = startcol; j < endcol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d5 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d5 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }
  *converge = (d1 + d2 + d3 + d4 + d5);  
  MPI_Waitall(8, req, stat);

//if (rank == 0){
  for(i = 1; i < lrows-1 ; i++)
  {
    for(j = 1; j < lcols-1 ; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }//}
  PrintGrid(mat1,nrows,ncols);

  //swap mat1 = mat2
/*  tmp = mat1;
  mat1 = mat2;
  mat2 = tmp;
*/
}


int main(int argc, char *argv[])
{

  int i, j;
  double **mat1, **mat2;
  //size of grid
  int nrows = 12, ncols = 12;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 10;
  int u, d, l, r;
  clock_t start, end;

  P = 2;
  Q = 2;

  lrows = 2 + nrows/Q;
  lcols = 2 + ncols/P;

  

  //MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .001;
  //dt = 5;
  mat1 = calloc(lrows, sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat1[i] = calloc(lcols, sizeof(double));
	

  mat2 = calloc(lrows, sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat2[i] = calloc(lcols, sizeof(double));

  InitGrid(mat1);
  InitGrid(mat2);
  time = 0;

  start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
      time = time + dt;
      Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy);
      iter++;
      if (converge < epsilon)
       	break;
  }
  end = clock();
  if (rank == 0)
  {
    printf("num of iterations: %d, with change of %lf\n", iter, converge);
    printf("Total time: %f seconds\n", (((double) end) - start)/CLOCKS_PER_SEC);
  }
  MPI_Finalize();

}
