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
#define EPSILON 0.001 /* threshold for convergence */

struct thread_args {
  long id;
  int i_init;
  int i_fin;
  int j_init;
  int j_fin;
};

void *Worker(void *);
void InitializeGrids();
void Barrier();

struct tms buffer;        /* used for timing */
clock_t start, finish;

pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived = 0;       /* count of the number who have arrived */

int gridSize, numWorkers, numIters, stripSize;
double maxDiff[MAXWORKERS];
double grid1[MAXGRID][MAXGRID];
double grid2[MAXGRID][MAXGRID];
double maxdiff;

/* main() -- read command line, initialize grids, and create threads
             when the threads are done, print the results */

int main(int argc, char* argv[]) {
  /* thread ids and attributes */
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;
  int i, j;
  //long li;
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
  
  stripSize = gridSize/numWorkers;
  InitializeGrids();

  unsigned int quadrants = 1;
  unsigned int quadSize = gridSize;

  while(!(quadSize & 1) && (quadrants << 2) <= numWorkers){
    quadrants <<= 2;
    quadSize >>= 1;
    //printf("quadrant gridsize = %d, num of quadrants = %d\n", quadSize, quadrants);
  }

  /*
  int x = quadrants - numWorkers;
  if(x >= 0)
    printf("leftover quadrants: %d\n", x);
  else{
    printf("leftover threads: %d\n", -x);
    //numWorkers += x;
  }
  */
  
  printf("Unused threads: %d\n", numWorkers - quadrants);
  numWorkers += quadrants - numWorkers;

  /* get indices of array */

  struct thread_args t_args[numWorkers];

  start = times(&buffer);
  /* create the workers, then wait for them to finish */
  printf("Create the worker threads\n");
  for (i = 0; i < numWorkers; i++)
    {
      //li = i;
      t_args[i].id = i;
      t_args[i].i_init = i * quadSize + 1;
      t_args[i].i_fin = t_args[i].i_init + quadSize - 1;
      t_args[i].j_init = i * quadSize + 1;
      t_args[i].j_fin = t_args[i].i_init + quadSize - 1;
      
      pthread_create(&workerid[i], &attr, Worker, (void *) &(t_args[i]));
    }
  printf("Wait for all worker threads to finish\n");
  for (i = 0; i < numWorkers; i++)
    pthread_join(workerid[i], NULL);

  finish = times(&buffer);
  /* print the results */
  
  printf("number of iterations:  %d\nmaximum difference:  %e\n",
          numIters, maxdiff);
  printf("start:  %ld   finish:  %ld\n", start, finish);
  printf("elapsed time:  %ld\n", finish-start);
  results = fopen("results", "w");
  for (i = 1; i <= gridSize; i++) {
    for (j = 1; j <= gridSize; j++) {
      fprintf(results, "%f ", grid2[i][j]);
    }
    fprintf(results, "\n");
  }
  
}


/* Each Worker computes values in one strip of the grids.
   The main worker loop does two computations to avoid copying from
   one grid to the other.  */

void *Worker(void *arg) {
  struct thread_args* t_arg = arg;
  long myid =(long) t_arg->id;
  double localDiff, temp;
  int i, j;
  //int first, last;
  int i_last = t_arg->i_fin;
  int j_last = t_arg->j_fin;

  printf("worker %ld (pthread id %lu has started\n", myid, pthread_self());

  /* determine first and last rows of my strip of the grids */
  //first = myid * stripSize + 1;
  //last = first + stripSize - 1;
  
  while (maxdiff > EPSILON) {
    /* update my points */
    for (i = t_arg->i_init; i <= i_last; i++) {
      for (j = t_arg->j_init; j <= j_last; j++) {
        grid2[i][j] = (grid1[i-1][j] + grid1[i+1][j] + 
          grid1[i][j-1] + grid1[i][j+1]) * 0.25;
      }
    }
    Barrier();
    /* update my points again */
    localDiff = 0.0;
    for (i = t_arg->i_init; i <= i_last; i++) {
      for (j = t_arg->j_init; j <= j_last; j++) {
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
