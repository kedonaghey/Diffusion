#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void InitFineGrid(double** mat, int nrows, int ncols)
{
  int i;
  int top = 20;
  int left = 50;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows; i+=2)
    {
      mat[i][0] = top;
    }
  for (i = 1; i < nrows; i+=2)
    {
      mat[i][0] = top/2;
    }
  //left and right rows initial conditions
  for ( i = 0; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left;//100
    }
  for (i=1; i < ncols; i+=2)
    {
      mat[nrows-1][i] = left/2;
    }

}

void InitGrid(double** mat, int nrows, int ncols, int crs, int ccs)
{
  int i;
  //bottom and top rows initial conditions
  for ( i = 0; i < nrows - crs +1; i++)
    {
      mat[i][0] = 20;//0
    }
  for (i = 0; i <nrows; i++)
      mat[i][ncols - 1] = 30;//100
  //left and right rows initial conditions
  for ( i = ccs -1; i < ncols; i++)
    {
      mat[nrows - 1][i] = 50;//0
    }
  for ( i = 0; i < ncols; i++)
      mat[0][i] = 40;//100

}

void Printgrid(double** mat, int nrows, int ncols)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }


}
void PrintGrid(double** mat, int nrows, int ncols, FILE *fp)
{
  int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%f ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }


}

void InjectCoarsePoints(double**mat, int nrows, int ncols, int crs, int ccs)
{
  int i, j;
  //top
  for (i = 1; i < ccs; i++)
      mat[nrows-crs][i] = 50;//0
  //right
  for(i = nrows-2; i > nrows - crs; i--)
      mat[i][ccs - 1] = 60;//100
  //corner
  mat[nrows-crs][ccs-1] = 100;
}

void InjectFinePoints(double** mat, int nrows, int ncols)
{
  int i, j;
  //top
  for (i = 1; i < ncols - 1; i++)
      mat[0][i] = 40;//0
  //right
  for(i = 1; i < nrows - 1; i++)
      mat[i][ncols - 1] = 30;//100
  //corner
  //double check
  mat[0][ncols-1] = 1;
}

/****    Interface Points *****/


double Compute_CornerTimeStep(double** mat1_refine, double** mat1, int i, int j, double* converge, double dt, double dx, double dy)
{

  double wxy = 0;
  double k = .02, diff = 0;
  wxy = k * dt/(dx*dx);

  mat2_refine[i][j] = /*prev point*/mat1_refine[i][j] + wxy * (16/15) * (/*prev point*/-4*mat1_refine[i][j] + .5 * /*bound*/mat1_refine[i][j-2]/
  + /*course*/mat1_refine[i+2][j] + .5 * /*course*/mat1_refine[i-2][j] +/*bound*/ mat1_refine[i][j+4] +/*mesh*/ mat1_refine[i-1][j-1]);

 return mat2_refine[i][j];
}


void Compute_InterfaceTimeStep(double** mat1_refine, double** mat1, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{

  int i,j;
  double wxy = 0;
  double k = .02, diff = 0;
  wxy = k * dt/(dx*dx);

  for (i = 1; i < nrows - 1; i++)
  {
    for (j = 1; j < ncols - 1; j++)
    {
      mat2_refine[i][j] = mat1_refine[i][j] + mat1_refine[i][j] * wxy * (dx/2) + wxy * (4/3) * (- 4*mat1_refine[i][j] + mat1_refine[i+2][j]  + .5 * mat1_refine[i][j+2]/
      + .5 * mat1_refine[i][j-2] + .5 * mat1_refine[i-1][j+1] + mat1_refine[i-1][j] + .5 * mat1_refine[i-1][j+1]);
    }
  }


}




void Compute_FineTimeStep(double** mat1_refine, double** mat2_refine, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{

  int i,j;
  double wxy = 0;
  double k = .02, diff = 0;
  int refinement = 2;
  double sxy = (double)dx/refinement;
  wxy = k * dt/(dx*dx);

  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
      mat2_refine[i][j] = mat1_refine[i][j] + wxy * sxy * sxy * mat1_refine[i][j] + wxy * ( .5 * (mat1_refine[i+1][j+1] + mat1_refine[i+1][j-1])\
      -2 * mat1_refine[i][j] + mat1_refine[i-1][j]);
      if ( diff < fabs(mat2_refine[i][j] - mat1_refine[i][j]))
      {
        diff = fabs(mat1_refine[i][j] - mat2_refine[i][j]);
      }
    }
  }

  *converge = diff;  

  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
    {
      mat1_refine[i][j] = mat2_refine[i][j];
    }
  }

}

void Compute_TimeStep(double** mat1, double** mat2, int nrows, int ncols, double* converge, double dt, double dx, double dy)
{
  int i, j;
  double **tmp;
  double wx = 0, wy = 0;
  //diffusitivity coefficient
  double k = .02;

//printf("dx = %lf, dy = %lf, dt = %lf\n", dx, dy, dt);
  //weight in x direction
  wx = k * dt/(dx*dx);
  //weight in y direction
  wy = k * dt/(dy*dy);

  //stability
  // wx + wy must be less than .5
  if ( wx + wy > .5)
  {
    printf("Does not reach requirements for stability\n");
    printf("wx = %lf, wy = %lf\n", wx, wy);
    exit(1);
  }

  
  double diff = 0;
  for (i = 1; i < nrows -1; i++)
  {
    for (j = 1; j < ncols -1; j++)
    {
      mat2[i][j] = mat1[i][j] + wx * (mat1[i+1][j] - 2*mat1[i][j] + mat1[i-1][j]) + wy * (mat1[i][j+1] - 2 * mat1[i][j] + mat1[i][j-1]);
      if ( diff < fabs(mat2[i][j] - mat1[i][j]))
      {
        diff = fabs(mat1[i][j] - mat2[i][j]);
      }
    }
  }

  	*converge = diff;  

  for(i = 0; i < nrows; i++)
  {
    for(j = 0; j < ncols; j++)
    {
      mat1[i][j] = mat2[i][j];
    }
  }

}


int main(int argc, char *argv[])
{

  int i, j;
  double **mat1, **mat2;
  double **mat1_refine, **mat2_refine;
  //size of grid
  int nrows = 10, ncols = 10;
  //2x2 for mesh in a 4x4 array
  double refrows = 5, refcols = 5;
  double dx = 0, dy= 0, dt, time;
  double converge = 0, epsilon;
  int iter, max_iter = 100;
  clock_t start, end;
  int prnt = 1;
  char output_file[80] = "k.txt";
  FILE *fp;

  //width of space steps in x and y direction
  dx = 1 /(double)(nrows - 1);
  dy = 1 /(double)(ncols - 1);

  epsilon = .001;
  dt = .001;
  //dt = 5;
  mat1 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat1[i] = calloc(ncols, sizeof(double));
	

  mat2 = calloc(nrows, sizeof(double *));
  for (i = 0; i < nrows; i++)
    mat2[i] = calloc(ncols, sizeof(double));

  int crs = refrows/2 + 1;
  int ccs = refcols/2 + 1;
  InitGrid(mat1, nrows, ncols, crs, ccs);
  InitGrid(mat2, nrows, ncols, crs, ccs);

  mat1_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat1_refine[i] = calloc(refcols, sizeof(double));
	

  mat2_refine = calloc(refrows, sizeof(double *));
  for (i = 0; i < refrows; i++)
    mat2_refine[i] = calloc(refcols, sizeof(double));

  InitFineGrid(mat1_refine, refrows, refcols);
  InitFineGrid(mat2_refine, refrows, refcols); 

  time = 0;

  fp = fopen ( output_file, "w" );
  start = clock();
  for(iter = 0; iter < max_iter; iter++)
  {
    //move to seperate func?
    //while ( converge >= epsilon)
    //{
      time = time + dt;
      Compute_TimeStep(mat1, mat2, nrows, ncols, &converge, dt, dx, dy);
      Compute_FineTimeStep(mat1_refine, mat2_refine, refrows, refcols, &converge, dt, dx, dy);
      Compute_InterfaceTimeStep();
      Compute_CornerTimeStep();
      InjectCoarsePoints(mat1, nrows, ncols, crs, ccs);
      InjectFinePoints(mat1_refine, refrows, refcols);
      /*if (iter%100 == 0)
      { 
        printf("%d %f\n", iter, converge);
        fprintf(fp, "\n");
	    //prnt = 2 * prnt;
      }*/
      iter++;
      if (converge < epsilon)
       	break;
    //}
  }
  end = clock();

  Printgrid(mat1_refine,refrows,refcols);
  //fprintf ( fp, "%d\n", nrows );
  //fprintf ( fp, "%d\n", ncols );

  PrintGrid(mat1, nrows, ncols, fp);
  //fprintf(fp, "num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);
  //fprintf(fp, "Total time: %f seconds\n", (((double) end ) - start)/CLOCKS_PER_SEC);

  fclose(fp);
  printf("num of iterations: %d, with change of %lf %d\n", iter, converge, prnt);


}

