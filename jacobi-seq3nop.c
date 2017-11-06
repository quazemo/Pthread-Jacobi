/* Jacobi iteration plain vanilla sequential version

   usage on Linux:
     gcc jacobi-seq3nop.c  -o jacobi-seq
     jacobi x y epsilon
            columns rows

*/
#define IDX(x, i, j) ((i)*(x)+(j))
#define ULONG_MAX    18446744073709551615UL
#define num_runs 1000
#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include <limits.h>
#define FAILURE 0
#define SUCCESS 1
struct tms buffer, start_struc, finish_struc;        /* used for timing */
clock_t start, finish, time_sample, least_time;
unsigned long least_run ;
int x, y, xm1, ym1;
double *grid1;
double *grid2;
double *startgrid1;
double *startgrid2;
double maxdiff; 
double eps;
int numIters=0;
unsigned long cycles, start_tsc, end_tsc;
unsigned long cycles_table[num_runs];
clock_t time_table[num_runs];
static __inline__ unsigned long Gettsc(void)
{
  unsigned a, d;
  asm ("cpuid");
  asm volatile("rdtsc" : "=a" (a), "=d" (d)); 
  return ((unsigned long)a) | (((unsigned long)d) << 32); 
}

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
  //  printf("Worker\n");
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
    // printf("maxdiff=%f, Iterations=%d\n",maxdiff,numIters);
  }
}

void ReInitializeGrids(double *grid1, double *startgrid1, double *grid2, double *startgrid2)
{
  /* re initialize grid1 and grid2 */
  
  int i, j;
  for (i = 0; i < y; i++)
    for (j = 0; j < x; j++) {
      grid1[IDX(x,i,j)]=startgrid1[IDX(x,i,j)];
      grid2[IDX(x,i,j)]=startgrid2[IDX(x,i,j)];
    }
}


void InitializeGrids(double *grid1, double *startgrid1, double *grid2, double *startgrid2)
{
  int i, j;
  double d;
  for (i = 0; i < y; i++)
    for (j = 0; j < x; j++) {
      scanf("%lg",&d);
      grid1[IDX(x,i,j)]=d;
      startgrid1[IDX(x,i,j)]=d;
      /* put the boundary values into the second grid as well 
	 they never change      */
      if ((j==0) || (j==xm1) || (i==0) || (i==ym1))
	{
	  grid2[IDX(x,i,j)]=grid1[IDX(x,i,j)];
	  startgrid2[IDX(x,i,j)]=grid1[IDX(x,i,j)];
	}
      else
	{
	  grid2[IDX(x,i,j)]=0.0;
	  startgrid2[IDX(x,i,j)]=0.0;
	}                                    // don't really need to do this but I want clean numbers on
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
  int runs, i, j, gridsize;
  int badParm=0;
  double *grid1;
  double *grid2;
  /* read command line and initialize grids */
  if (argc < 3 ) {
    printf("Need three parameters: x y epsilon\n");
    printf("Bye\n");
    exit(FAILURE);
  }
  x = atoi(argv[1]);
  y = atoi(argv[2]);
  eps = atof(argv[3]);
  if (x <=0)
    { printf("x needs to be greater than 0");
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
 startgrid1=(double *) malloc(gridsize);
 startgrid2=(double *) malloc(gridsize);
 
 InitializeGrids(grid1,startgrid1,grid2,startgrid2);
 // PrintGrid(grid1);
 // printf("\n");
 // PrintGrid(grid2);
 for (runs=0; runs<num_runs; runs++)
   {
     maxdiff=99999.9;
     times(&start_struc);
     start_tsc=Gettsc();
     Worker(grid1,grid2);
     end_tsc=Gettsc();
     times(&finish_struc);
     /* print the results */
     //     printf("number of iterations:  %d\nmaximum difference:  %e\n",numIters, maxdiff);
     // printf("start:  %ld   finish:  %ld\n", start, finish);
     // printf("elapsed time:  %ld\n", finish-start);
     // need code to handle wrap around
     // printf("start tsc:  %ld   end tsc:  %ld\n", start_tsc, end_tsc);
     cycles= (end_tsc > start_tsc)? end_tsc-start_tsc : (ULONG_MAX-start_tsc + end_tsc) ;
     cycles_table[runs]=cycles;
     time_sample=finish_struc.tms_utime - start_struc.tms_utime;
     time_table[runs]=time_sample;
     ReInitializeGrids(grid1,startgrid1,grid2,startgrid2);
     // printf("number of cycles:  %ld\n", cycles);
   }
 least_run=ULONG_MAX;
 least_time=99999;
 for (runs=0; runs<num_runs; runs++)
   {
     printf("number of cycles:  %ld\n", cycles_table[runs]);
     printf("utime:  %ld\n", time_table[runs]);
     least_run=(cycles_table[runs] < least_run) ? cycles_table[runs] : least_run;
     least_time=(time_table[runs] < least_time) ? time_table[runs] : least_time;
   }
 printf("lowest number of cycles=%ld\n",least_run);
 printf("lowest time=%ld\n",least_time);
 
 // PrintGrid(grid1);
 exit(SUCCESS);
}

