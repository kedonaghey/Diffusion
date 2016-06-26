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
    for(i=0;i<lrows;i++)
	mat[i][0] = 20;
  }
  if(rank%P == P-1)
  {
    for(i=0;i<lrows;i++)
	mat[i][lcols-1] = 30;
  }
  if(rank/P == 0)
  {
    for(i=0;i<lcols;i++)
	mat[0][i] = 50;
  }
  if(rank/P == Q-1)
  {
    for(i=0;i<lcols;i++)
	mat[lrows - 1][i] = 40;
  }
/*
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i++)
    {
      mat[i][0] = 20;//0
      mat[i][ncols - 1] = 30;//100
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i++)
    {
      mat[0][i] = 50;//100
      mat[nrows - 1][i] = 40;//0
    }
*/
}
void PrintGrid(double** mat, int nrows, int ncols)
{
/*
 MPI_Datatype BLOCK;
 double *fullgrid;
 int i, j, k, l;
 int offset, x, y;
 if(rank == 0) {
  i = size*lrows*lcols*sizeof(double);
  fullgrid = malloc(i);
}
MPI_Type_vector(lrows, lcols, lcols, MPI_DOUBLE, &BLOCK);
	MPI_Type_commit(&BLOCK);

	// Easiest to collect the whole grid onto rank 0 and then print 

	MPI_Gather(&mat[0][0], 1, BLOCK, fullgrid, 1, BLOCK, 0, MPI_COMM_WORLD);

	// Now detangle the mess that is the full grid array 

	if(rank == 0) {

	// For each row of MPI tasks
		for(i=0;i<Q;i++) {
			// For each of the rows 
			for(j=0;j<lrows;j++) {
				// For each of the subblocks 
			for(k=0;k<P;k++) {
				//Print each element 
				for(l=0;l<lcols;l++) {
						offset = i*lrows*lcols*P + j*lcols + k*lrows*lcols + l;	
						//printf("[%d]\t", offset);
						printf("%.1f\t", fullgrid[offset]);
						}
						printf("|");
						
						}
						printf("\n");
						}
						}
						//free(fullgrid);
						}
			//MPI_Type_free(&BLOCK);
*/
int i, j, k, l;
  for(l=0;l<Q;l++) {
    for (i = 0; i < lrows; i++) {
	for(k=0;k<P;k++) {
        for (j = 0; j < lrows; j++) {
            printf("%f ", mat[i][j]);
        }
        }
        printf("\n");
      }
}
    

}
void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{
  int i, j;
  double **tmp;
  double wx = 0, wy = 0;
  //diffusitivity coefficient - can equal .002
  double k = .02;
  double d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0, d5 = 0.0;

//printf("dx = %lf, dy = %lf, dt = %lf\n", dx, dy, dt);
  //weight in x direction
  wx = k * dt/(dx*dx);
  //weight in y direction
  wy = k * dt/(dy*dy);

//printf("wx = %lf, wy = %lf\n", wx, wy);
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
  int srow, scol, erow, ecol;

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



/*  if (rank%P == 0)
    scol = 2;
  else
    scol = 1;
  if (rank%P == P - 1)
    ecol = lcols - 2;
  else
    ecol = lcols - 1;
  if (rank/P == 0)
    srow = 1;
  else
    srow = 1;
  if (rank/P == Q - 1)
    erow = lrows - 2;
  else
    erow = lrows - 1;
*/
  if (rank%P == 0)
    scol = 1;
  else
    scol = 0;
  if (rank%P == P - 1)
    ecol = lcols -1;
  else
    ecol = lcols ;
  if (rank/P == 0)
    srow = 1;
  else
    srow = 0;
  if (rank/P == Q - 1)
    erow = lrows - 1;
  else
    erow = lrows ;
  //solving edges
  if (rank/P != 0)
  {
    //not in top blocks
    i = 0;
    for(j = scol; j < ecol; j++)
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
    i = lrows - 1;
    for(j = scol; j < ecol; j++)
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
    j = 0;
    for(i = srow; i < erow; i++)
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
    j = lcols - 1;
    for(i = srow; i < erow; i++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d4 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d4 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  for(i = srow; i < erow; i++)
  {
    for(j = scol; j < ecol; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( d5 < fabs(mat2[i][j] - mat1[i][j]))
      {
        d5 = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }
  *converge = (d1 + d2 + d3 + d4 + d5);  
  //halo exchange
  // Send rows up and down
  MPI_Irecv(mat2[0], lcols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[0]);
  MPI_Irecv(mat2[lrows-1], lcols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[1]);
  MPI_Isend(mat2[lrows-2], lcols, MPI_DOUBLE, d, 0, MPI_COMM_WORLD, &req[2]);
  MPI_Isend(mat2[1], lcols, MPI_DOUBLE, u, 0, MPI_COMM_WORLD, &req[3]);
  // Send columns left and right
  MPI_Irecv(&mat2[0][0], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[4]);
  MPI_Irecv(&mat2[0][lcols-1], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[5]);
  MPI_Isend(&mat2[0][1], 1, MPI_CLM, l, 0, MPI_COMM_WORLD, &req[6]);
  MPI_Isend(&mat2[0][lcols-2], 1, MPI_CLM, r, 0, MPI_COMM_WORLD, &req[7]);

/*
  double diff = 0;
  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( diff < fabs(mat2[i][j] - mat1[i][j]))
      {
        diff = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }
*/
  //solving inner

  MPI_Waitall(8, req, stat);


/*  
  double diff = 0;
  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( diff < fabs(mat2[i][j] - mat1[i][j]))
      {
        diff = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }*/

  for(i = 0; i < lrows; i++)
  {
    for(j = 0; j < lcols; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }

 // PrintGrid(mat1,nrows,ncols);
/*  tmp = mat1;
  mat1 = mat2;
  mat2 = tmp;*/
}


int main(int argc, char *argv[])
{

  int i, j;
  double **mat1, **mat2;
  //size of grid
  int nrows = 10, ncols = 10;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 10;
  int u, d, l, r;
  clock_t start, end;

  P = 1;
  Q = 1;

  lrows =/* 2 +*/ nrows/Q;
  lcols =/* 2 +*/ ncols/P;

  

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
/*  mat1 = malloc(lrows * sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat1[i] = malloc(lcols * sizeof(double));
	

  mat2 = malloc(lrows * sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat2[i] = malloc(lcols * sizeof(double));

  InitGrid(mat1);
  InitGrid(mat2);
*/

/*
  mat1 = malloc(nrows *sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = malloc(ncols *sizeof(double));
	

  mat2 = malloc(nrows *sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = malloc(ncols *sizeof(double));
*/

  mat1 = calloc(lrows, sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat1[i] = calloc(lcols, sizeof(double));
	

  mat2 = calloc(lrows, sizeof(double *));
  for (i = 0; i < lrows; i++)
    mat2[i] = calloc(lcols, sizeof(double));

  InitGrid(mat1);
  InitGrid(mat2);
  //if (rank == 0)
  //PrintGrid(mat1,nrows,ncols);
  time = 0;

  start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
    //move to seperate func?
    //while ( converge >= epsilon)
    //{
      time = time + dt;
      Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy);
      iter++;
      if (converge < epsilon)
       	break;
    //}
  }
  end = clock();
  if (rank == 0)
  {
    printf("num of iterations: %d, with change of %lf\n", iter, converge);
    printf("Total time: %f seconds\n", (((double) end) - start)/CLOCKS_PER_SEC);
  }
  MPI_Finalize();

}
