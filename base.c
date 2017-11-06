/* Jacobi iteration using pthreads

   usage on Linux:
     gcc -g base.c -lpthread -o base
     base gridSize numWorkers
     
     Brendan Baalke
     Iris Larsen
     Alexander Lee
     CSCI 415, Fall Quarter, 2017
  
     Base code provided by Dr. Michael Meehan;
     edited for experimentation purposes

     utilizes EPSILON instead of numIters
*/

#define _REENTRANT
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <sys/times.h>
#include <limits.h>
#define SHARED 1
#define MAXGRID 7000   /* maximum grid size, including boundaries */
#define MAXWORKERS 20  /* maximum number of worker threads */
#define EPSILON 0.0001 /* threshold for convergence */
#define num_runs 1000

static __inline__ unsigned long Gettsc(void);
void *Worker(void *);
void InitializeGrids();
void Barrier();

struct tms buffer, start_struc, finish_struc;        /* used for timing */
clock_t start, finish, time_sample, least_time;
unsigned long least_run ;

unsigned long cycles, start_tsc, end_tsc;
unsigned long cycles_table[num_runs];
clock_t time_table[num_runs];

pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived = 0;       /* count of the number who have arrived */

int gridSize, numWorkers, numIters, stripSize;
double maxdiff;
double maxDiff[MAXWORKERS];
double grid1[MAXGRID][MAXGRID], grid2[MAXGRID][MAXGRID];


/* main() -- read command line, initialize grids, and create threads
             when the threads are done, print the results */

int main(int argc, char *argv[]) {
  /* thread ids and attributes */
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;
  int i, j, runs;
  long li;
  numIters = 0;
  maxdiff = EPSILON + 1;
  FILE *results;

  /* set global thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  /* initialize mutex and condition variable */
  pthread_mutex_init(&barrier, NULL);
  pthread_cond_init(&go, NULL);

  /* read command line and initialize grids */
  gridSize = atoi(argv[1]);
  numWorkers = atoi(argv[2]);

  printf("gridSize: %d, number of threads: %d\n", gridSize, numWorkers);
  
  stripSize = gridSize / numWorkers;
  InitializeGrids();

  for (runs=0; runs<num_runs; runs++)
    {
      maxdiff=99999.9;
      times(&start_struc);
      start_tsc=Gettsc();
  
      /* create the workers, then wait for them to finish */
      //printf("Create the worker threads\n");
      for (i = 0; i < numWorkers; i++) {
	li=i;
	pthread_create(&workerid[i], &attr, Worker, (void *) li);
      }
      //printf("Wait for all worker threads to finish\n");
      for (i = 0; i < numWorkers; i++)
	pthread_join(workerid[i], NULL);

      end_tsc=Gettsc();
      times(&finish_struc);

      cycles= (end_tsc > start_tsc)? end_tsc-start_tsc : (ULONG_MAX-start_tsc + end_tsc) ;
      cycles_table[runs]=cycles;
      time_sample=finish_struc.tms_utime - start_struc.tms_utime;
      time_table[runs]=time_sample;

      InitializeGrids();
    }

  /* print the results 
  for (i = 0; i < numWorkers; i++)
    if (maxdiff < maxDiff[i])
      maxdiff = maxDiff[i];
  printf("number of iterations:  %d\nmaximum difference:  %e\n",
          numIters, maxdiff);
  */

  least_run=ULONG_MAX;
  least_time=99999;
  for (runs=0; runs<num_runs; runs++)
    {
      //printf("number of cycles:  %lu\n", cycles_table[runs]);
      //printf("utime:  %ld\n", time_table[runs]);
      least_run=(cycles_table[runs] < least_run) ? cycles_table[runs] : least_run;
      least_time=(time_table[runs] < least_time) ? time_table[runs] : least_time;
    }
  printf("lowest number of cycles=%lu\n",least_run);
  printf("lowest time=%ld\n",least_time);
  

  /*
  results = fopen("results", "w");
  for (i = 1; i <= gridSize; i++) {
    for (j = 1; j <= gridSize; j++) {
      fprintf(results, "%f ", grid2[i][j]);
    }
    fprintf(results, "\n");
  }
  */
}


/* Each Worker computes values in one strip of the grids.
   The main worker loop does two computations to avoid copying from
   one grid to the other.  */

void *Worker(void *arg) {
  long myid =(long)arg;
  double localDiff, temp;
  int i, j, iters;
  int first, last;

  //printf("worker %ld (pthread id %lu has started\n", myid, pthread_self());

  /* determine first and last rows of my strip of the grids */
  first = myid * stripSize + 1;
  //last = first + stripSize - 1;

  if(myid == numWorkers-1)
    last = gridSize;
  else
    last = first + stripSize - 1;
  
  while (maxdiff > EPSILON) {
    /* update my points */
    for (i = first; i <= last; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid2[i][j] = (grid1[i-1][j] + grid1[i+1][j] + 
          grid1[i][j-1] + grid1[i][j+1]) * 0.25;
      }
    }
    Barrier();
    /* update my points again */
    localDiff = 0.0;
    for (i = first; i <= last; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid1[i][j] = (grid2[i-1][j] + grid2[i+1][j] +
          grid2[i][j-1] + grid2[i][j+1]) * 0.25;
        temp = grid1[i][j] - grid2[i][j];
        if (temp < 0) {
          temp *= -1;
        }
        if (localDiff < temp) {
          localDiff = temp;
        }
      }
    }
    /* compute the maximum difference in my strip and set global variable */
    maxDiff[myid] = localDiff;
    Barrier();
    if (myid == 0) {
      maxdiff = 0.0;
      numIters++;
      for (i = 0; i < numWorkers; i++) {
        if (maxdiff < maxDiff[i]) {
          maxdiff = maxDiff[i];
        }
      }
    }
    Barrier();
  }
  return NULL;
}

void InitializeGrids() {
  /* initialize the grids (grid1 and grid2)
     set boundaries to 1.0 and interior points to 0.0  */
  int i, j;
  for (i = 0; i <= gridSize+1; i++)
    for (j = 0; j <= gridSize+1; j++) {
      grid1[i][j] = 0.0;
      grid2[i][j] = 0.0;
    }
  for (i = 0; i <= gridSize+1; i++) {
    grid1[i][0] = 1.0;
    grid1[i][gridSize+1] = 1.0;
    grid2[i][0] = 1.0;
    grid2[i][gridSize+1] = 1.0;
  }
  for (j = 0; j <= gridSize+1; j++) {
    grid1[0][j] = 1.0;
    grid2[0][j] = 1.0;
    grid1[gridSize+1][j] = 1.0;
    grid2[gridSize+1][j] = 1.0;
  }
}

void Barrier() {
  pthread_mutex_lock(&barrier);
  numArrived++;
  if (numArrived == numWorkers) {
    numArrived = 0;
    pthread_cond_broadcast(&go);
  } else
    pthread_cond_wait(&go, &barrier);
  pthread_mutex_unlock(&barrier);
}

static __inline__ unsigned long Gettsc(void)
{
  unsigned a, d;
  asm ("cpuid");
  asm volatile("rdtscp" : "=a" (a), "=d" (d)); 
  return ((unsigned long)a) | (((unsigned long)d) << 32); 
}
