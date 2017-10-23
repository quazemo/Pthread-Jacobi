/* Jacobi iteration plain vanilla sequential version

   usage on Linux:
     gcc -g -Wall jacobi-seq.c -o jacobi-seq
     jacobi length width epsilon

   Brendan Baalke
   Iris Larsen
   Alexander Lee
   CSCI 415, Fall Quarter, 2017

   Base code provided by Dr. Michael Meehan;
   edited for experimentation purposes

   Modified 10/23/2017 by Alexander Lee

*/
#define IDX(x, i, j) ((i)*(x)+(j))

#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include <limits.h>
#include <time.h>
#define FAILURE 0
#define SUCCESS 1
struct tms buffer;        /* used for timing */
//clock_t start, finish;

struct timespec start, finish;

int x, y, xm1, ym1;
double *grid1;
double *grid2;
double maxdiff = 99999.9;
double eps;
int numIters=0;

void PrintGrid(double *grid)
{
  int i,j;
  //  printf("PrintGrid\n");
  for (i = 0; i < y; i++) {
    for (j = 0; j < x; j++) {
      printf("%f ", grid[IDX(x,i,j)]);
    }
    printf("\n");
  }
}

void Worker(double *grid1, double *grid2 )
{
  double temp;
  int i, j;
  printf("Worker\n");
  while (maxdiff > eps) {
    for (i = 1; i < ym1; i++) {
      for (j = 1; j < xm1; j++) {
        grid2[IDX(x,i,j)] = (grid1[IDX(x,i-1,j)] + grid1[IDX(x,i+1,j)] + 
			     grid1[IDX(x,i,j-1)] + grid1[IDX(x,i,j+1)]) * 0.25;
      }
    }
    numIters++;
    for (i = 1; i < ym1; i++) {
      for (j = 1; j < xm1; j++) {
        grid1[IDX(x,i,j)] = (grid2[IDX(x,i-1,j)] + grid2[IDX(x,i+1,j)] + 
			     grid2[IDX(x,i,j-1)] + grid2[IDX(x,i,j+1)]) * 0.25;
      }
    }
    numIters++;
  
  /* compute the maximum difference into global variable */
    maxdiff=0;
  for (i = 1; i < ym1; i++) {
    for (j = 1; j < xm1; j++) {
      temp = grid1[IDX(x,i,j)]-grid2[IDX(x,i,j)];
      if (temp < 0)
        temp = -temp;
      if (maxdiff <  temp)
        maxdiff = temp;
    }
  }
  //printf("maxdiff=%f, Iterations=%d\n",maxdiff,numIters);
  }
}

void InitializeGrids(double *grid1, double *grid2)
{
  /* initialize the grids (grid1 and grid2) */

  int i, j;
  double d = 1.0;
  //  printf("Initialize Grids\n");
  //scanf("%lg",&d);
  for (i = 0; i < y; i++)
    for (j = 0; j < x; j++) {
      //      printf("%lg\n",d);
      /* put the boundary values into the second grid as well 
	 they never change      */
      if ((j==0) || (j==xm1) || (i==0) || (i==ym1))
	{
	  grid1[IDX(x,i,j)]=d;
	  grid2[IDX(x,i,j)]=grid1[IDX(x,i,j)];
	}
      else
	{
	  grid2[IDX(x,i,j)]=0.0;
	} // don't really need to do this but I want clean numbers on
                             // debug output
    }
}

/* main() 
   read command line, 
   initialize grids, 
   run till convergence             
   print the results 
*/

int main(int argc, char *argv[])
{
  //int i, j, gridsize;
  int gridsize;
  int badParm=0;
  double *grid1;
  double *grid2;
  /* read command line and initialize grids */
  if (argc < 3 ) {
    printf("Need three parameters: x y epsilon\n");
    printf("Bye\n");
    exit(FAILURE);
  }
  x = atoi(argv[1]) + 2;
  y = atoi(argv[2]) + 2;
  eps = atof(argv[3]);
  if (x <=0)
    { printf("length needs to be greater than 0");
      badParm=1;
    }
  if (y <=0)
    { printf("y needs to be greater than 0");
      badParm=1;
    }
 if (eps <=0)
    { printf("eps needs to be greater than 0");
      badParm=1;
    }
 if (badParm) exit(FAILURE);
 printf("x=%d, y=%d, epsilon=%f\n",x,y,eps);
 xm1=x-1;
 ym1=y-1;
 gridsize=x*y*sizeof(double);
 grid1=(double *) malloc(gridsize);
 grid2=(double *) malloc(gridsize);
 InitializeGrids(grid1,grid2);
 //PrintGrid(grid1);
 //printf("\n");
 //PrintGrid(grid2);
 //start = times(&buffer);

 clock_gettime(CLOCK_MONOTONIC, &start);

 Worker(grid1,grid2);

 clock_gettime(CLOCK_MONOTONIC, &finish);
 //finish = times(&buffer);
 /* print the results */
 printf("number of iterations:  %d\nmaximum difference:  %e\n",numIters, maxdiff);
 //printf("start:  %ld   finish:  %ld\n", start, finish);
 //printf("elapsed time:  %ld\n", finish-start);
  printf("start:  %f   finish:  %f\n", start.tv_sec + (start.tv_nsec/1000000000.0), finish.tv_sec + (finish.tv_nsec/1000000000.0));
  printf("elapsed time:  %f seconds\n", (double) (finish.tv_sec-start.tv_sec) + ((finish.tv_nsec - start.tv_nsec)/1000000000.0));
 //PrintGrid(grid1);
 exit(SUCCESS);
}

