/* Jacobi iteration using pthreads

   usage on Linux:
     gcc jacobi.c -lpthread -o jacobi
     jacobi gridSize numWorkers
  
  Brendan Baalke
  Iris Larsen
  Alexander Lee
  CSCI 415, Fall Quarter, 2017
  
  Base code provided by Dr. Michael Meehan;
  edited for experimentation purposes
  
  This version directly allocates 2 blocks of
  memory from the heap, then uses some arithmetic
  to access data in a single pointer dereference instead
  of 2 dereferences.
*/

#define _REENTRANT
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <sys/times.h>
#include <limits.h>
#define SHARED 1
#define MAXGRID 10000   /* maximum grid size, including boundaries */
#define MAXWORKERS 20  /* maximum number of worker threads */
#define EPSILON 0.00001 /* threshold for convergence */

void *Worker(void *);
void InitializeGrids();
void Barrier();

struct tms buffer;        /* used for timing */
//clock_t start, finish;
struct timespec start, finish;

pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived = 0;       /* count of the number who have arrived */

int gridSize, numWorkers, numIters, stripSize;
double maxDiff[MAXWORKERS];
double* grid1;
double* grid2;
double maxdiff;


/* main() -- read command line, initialize grids, and create threads
             when the threads are done, print the results */

int main(int argc, char* argv[]) {
  /* thread ids and attributes */
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;
  int i, j;
  long li;
  maxdiff = EPSILON + 1;
  numIters = 0;
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
  
  grid1 = malloc(sizeof(double) * (gridSize + 2) * (gridSize + 2));
  grid2 = malloc(sizeof(double) * (gridSize + 2) * (gridSize + 2));
  
  stripSize = gridSize / numWorkers;
  InitializeGrids();

  //start = times(&buffer);
  clock_gettime(CLOCK_MONOTONIC, &start);
  
  printf("gridSize: %d, number of threads: %d\n", gridSize, numWorkers);
  
  /* create the workers, then wait for them to finish */
  printf("Create the worker threads\n");
  for (i = 0; i < numWorkers; i++)
    {
      li = i;
      pthread_create(&workerid[i], &attr, Worker, (void *) li);
    }
  printf("Wait for all worker threads to finish\n");
  for (i = 0; i < numWorkers; i++)
    pthread_join(workerid[i], NULL);

  //finish = times(&buffer);
  clock_gettime(CLOCK_MONOTONIC, &finish);
  
  /* print the results */
  printf("number of iterations:  %d\nmaximum difference:  %e\n",
          numIters, maxdiff);
  //printf("start:  %ld   finish:  %ld\n", start, finish);
  //printf("elapsed time:  %ld\n", finish-start);

  printf("start:  %f   finish:  %f\n", start.tv_sec + (start.tv_nsec/1000000000.0), finish.tv_sec + (finish.tv_nsec/1000000000.0));
  printf("elapsed time:  %f seconds\n", (double) (finish.tv_sec-start.tv_sec) + ((finish.tv_nsec - start.tv_nsec)/1000000000.0));
  
  results = fopen("results", "w");
  for (i = gridSize + 3; i < (gridSize + 2) * (gridSize + 2); i += gridSize + 2) {
    for (j = 0; j < gridSize; j++) {
      fprintf(results, "%f ", grid2[i + j]);
    }
    fprintf(results, "\n");
  }
  free(grid1);
  free(grid2);
}


/* Each Worker computes values in one strip of the grids.
   The main worker loop does two computations to avoid copying from
   one grid to the other.  */

void *Worker(void *arg) {
  long myid = (long) arg;
  double localDiff, temp;
  int i, j;
  int first, last;

  printf("worker %ld (pthread id %lu has started\n", myid, pthread_self());

  /* determine first and last rows of my strip of the grids */
  first = (myid * stripSize * (gridSize + 2)) + gridSize + 3;
  last = first + (stripSize * (gridSize + 2));
  
  while (maxdiff > EPSILON) {
    /* update my points */
    for (i = first; i < last; i += gridSize + 2) {
      for (j = 0; j < gridSize; j++) {
        grid2[i + j] = (
          grid1[i + j - gridSize - 2] + //one row up
          grid1[i + j + gridSize + 2] + //one row down
          grid1[i + j-1] + //one column left
          grid1[i + j + 1]) //one column right
          * 0.25; //average the four points
      }
    }
    Barrier();
    /* update my points again, reverse grids */
    localDiff = 0.0;
    for (i = first; i < last; i += gridSize + 2) {
      for (j = 0; j < gridSize; j++) {
        grid1[i + j] = (
          grid2[i + j - gridSize - 2] + //one row up
          grid2[i + j + gridSize + 2] + //one row down
          grid2[i + j-1] + //one column left
          grid2[i + j + 1]) //one column right
          * 0.25; //average the four points
        temp = grid1[i + j] - grid2[i + j];
        if (temp < 0) {
          temp *= -1;
        }
        if (localDiff < temp) {
          localDiff = temp;
        }
      }
    }
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
  int i;
  for (i = 0; i < (gridSize + 2)  * (gridSize + 2); i++) {
      grid1[i] = 0.0;
      grid2[i] = 0.0;
  }
  for (i = 0; i < gridSize + 2; i++) {
    grid1[i] = 1.0; //top row
    grid1[((gridSize + 2) * (gridSize + 2)) - i - 1] = 1.0; //bottom row
    grid1[i * (gridSize + 2)] = 1.0; //far left column
    grid1[(i * (gridSize + 2)) + gridSize + 1] = 1.0; //far right column

    grid2[i] = 1.0;
    grid2[((gridSize + 2) * (gridSize + 2)) - i - 1] = 1.0;
    grid2[i * (gridSize + 2)] = 1.0;
    grid2[(i * (gridSize + 2)) + gridSize + 1] = 1.0;
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
