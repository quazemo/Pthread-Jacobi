/* Jacobi iteration using pthreads

   usage on Linux:
     gcc -g -Wall jacobi_box.c -lpthread -o jacobi_box
     jacobi gridSize numWorkers
  
  Brendan Baalke
  Iris Larsen
  Alexander Lee
  CSCI 415, Fall Quarter, 2017
  
  Base code provided by Dr. Michael Meehan;
  edited for experimentation purposes

  Modified 10/23/2017 by Alexander Lee
*/

#define _REENTRANT
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <sys/times.h>
#include <limits.h>
#include <string.h>
#define SHARED 1
#define MAXGRID 10000   /* maximum grid size, including boundaries */
#define MAXWORKERS 20  /* maximum number of worker threads */
#define EPSILON 0.001 /* threshold for convergence */

//Stolen from jacobi-deq2.c
#define IDX(x, i, j) ((i)*(x)+(j))

struct thread_args {
  long id;
  //int i_init;
  //int i_fin;
  int* j_init;
  int* j_fin;
};

void *Worker(void *);
void InitializeGrids();
void Barrier();

//struct tms buffer;        /* used for timing */
struct timespec start, finish;
//clock_t start, finish;

pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived = 0;       /* count of the number who have arrived */

int gridSize, numWorkers, numIters, stripSize;
double maxDiff[MAXWORKERS];

double* grid1;
double* grid2;

double maxdiff;

struct thread_args** t_args;


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
  gridSize = atoi(argv[1]) + 2;
  numWorkers = atoi(argv[2]);

  if(numWorkers > gridSize-2)
    numWorkers = gridSize-2;

  grid1 = malloc(sizeof(double*) * (gridSize) * (gridSize));
  grid2 = malloc(sizeof(double*) * (gridSize)  * (gridSize));
  
  stripSize = (gridSize-2)/numWorkers;
  InitializeGrids();

  //unsigned int quadrants = 1;
  //unsigned int quadSize = gridSize;

  /*
  while(!(quadSize & 1) && (quadrants << 2) <= numWorkers){
    quadrants <<= 2;
    quadSize >>= 1;
    //printf("quadrant gridsize = %d, num of quadrants = %d\n", quadSize, quadrants);
  }
  */
  //int i_i[];
  //int i_f[];
  //int j_i[numWorkers];
  //int j_f[numWorkers];

  //int boxes_per_row = stripSize;
  //int box_length = stripSize;
  int leftover = (gridSize-2)%numWorkers;
  
  t_args = (struct thread_args**) malloc(sizeof(struct thread_args*) * (numWorkers));
  for (int k = 0; k < numWorkers; k++){
    t_args[k] = (struct thread_args*) malloc(sizeof(struct thread_args) * (numWorkers));
    t_args[k]->j_init = malloc(sizeof(int) * (numWorkers+1));
    t_args[k]->j_fin = malloc(sizeof(int) * (numWorkers+1));
  }
  //t_args[numWorkers] = NULL;
  //1 extra for null terminated


  int k, j_i, j_f;
  
  for(k = 0, j_i = 1, j_f = 1; k < numWorkers; k++){
    if(k != 0)
      j_i = j_f+1;

    j_f = j_i+stripSize-1;
    //if(stripSize > 1)
    //  j_f += stripSize;
    
    if(leftover > 0){
      j_f++;
      leftover--;
    }
    t_args[0]->id = 0;
    t_args[0]->j_init[k] = j_i;
    t_args[0]->j_fin[k] = j_f;
    printf("k = %d, j_i = %d, j_f = %d\n", k, j_i, j_f);
  }
  for(k = 1; k < numWorkers; k++){
    t_args[k]->id = k;
    memcpy(&(t_args[k]->j_init[0]), &(t_args[0]->j_init[0]), sizeof(int) * numWorkers);
    memcpy(&(t_args[k]->j_fin[0]), &(t_args[0]->j_fin[0]), sizeof(int) * numWorkers);
  }

  for(k = 0; k < numWorkers; k++){
    printf("id = %ld\n", t_args[k]->id);
    for(int n = 0; n < numWorkers; n++)
      printf("j_i = %d j_f = %d\n", t_args[k]->j_init[n], t_args[k]->j_fin[n]);
    printf("\n");
  }

  clock_gettime(CLOCK_MONOTONIC, &start);
  /* create the workers, then wait for them to finish */
  printf("Create the worker threads\n");
  for (i = 0; i < numWorkers; i++)
    {
      //li = i;
      pthread_create(&workerid[i], &attr, Worker, (void *) t_args[i]);
    }
  printf("Wait for all worker threads to finish\n");
  for (i = 0; i < numWorkers; i++)
    pthread_join(workerid[i], NULL);

  clock_gettime(CLOCK_MONOTONIC, &finish);
  /* print the results */
  
  printf("number of iterations:  %d\nmaximum difference:  %e\n",
          numIters, maxdiff);
  printf("start:  %f   finish:  %f\n", start.tv_sec + (start.tv_nsec/1000000000.0), finish.tv_sec + (finish.tv_nsec/1000000000.0));
  printf("elapsed time:  %f seconds\n", (double) (finish.tv_sec-start.tv_sec) + ((finish.tv_nsec - start.tv_nsec)/1000000000.0));
  results = fopen("results", "w");
  for (i = 1; i < gridSize-1; i++) {
    for (j = 1; j < gridSize-1; j++) {
      fprintf(results, "%f ", grid2[IDX(gridSize,i,j)]);
    }
    fprintf(results, "\n");
  }

  for (int k = 0; k < numWorkers; k++){
    free(t_args[k]);
  }
  free(t_args);

  free(grid1);
  free(grid2);
  
}


/* Each Worker computes values in one strip of the grids.
   The main worker loop does two computations to avoid copying from
   one grid to the other.  */

void *Worker(void *arg) {
  struct thread_args* my_arg = arg; 
  long myid =(long) (my_arg->id);
  double localDiff, temp;
  int i, j, ix;
  int first, last;

  printf("worker %ld (pthread id %lu has started\n", myid, pthread_self());

  /* determine first and last rows of my strip of the grids 
  first = myid * stripSize + 1;

  if(myid == numWorkers-1)
    last = gridSize-1;
  else
    last = first + stripSize;
  */
  
  //first = my_arg->i_init;
  //last = my_arg->i_fin;
  
  while (maxdiff > EPSILON) {
    /* update my points */

    /* OLD
    for (i = first; i < last; i++) {
      for (j = 1; j < gridSize-1; j++) {
        grid2[IDX(gridSize,i,j)] = (grid1[IDX(gridSize,i-1,j)] + grid1[IDX(gridSize,i+1,j)] + 
				    grid1[IDX(gridSize,i,j-1)] + grid1[IDX(gridSize,i,j+1)]) * 0.25;
      }
    }
    */
    
    for (int k = 0; k < numWorkers; k++){
      first = my_arg->j_init[k];
      last = my_arg->j_fin[k];
      //printf("\n%d\n", k);

      ix = first*gridSize;
      //int fx = last*gridSize;
      for (i = first; i <= last; i++) {
	//int ix = i*gridSize;
	for (j = my_arg->j_init[myid]; j <= my_arg->j_fin[myid]; j++) {
	  grid2[ix+j] = (grid1[ix-gridSize+j] + grid1[ix+gridSize+j] + 
			 grid1[ix+j-1] + grid1[ix+j+1]) * 0.25;
	}
	ix+=gridSize;
      }
    }
    Barrier();
      /* update my points again */
    localDiff = 0.0;

      /* OLD
	 for (i = first; i < last; i++) {
	 for (j = 1; j < gridSize-1; j++) {
	 grid1[IDX(gridSize,i,j)] = (grid2[IDX(gridSize,i-1,j)] + grid2[IDX(gridSize,i+1,j)] + 
	 grid2[IDX(gridSize,i,j-1)] + grid2[IDX(gridSize,i,j+1)]) * 0.25;
	 temp = grid1[IDX(gridSize,i,j)] - grid2[IDX(gridSize,i,j)];
	 if (temp < 0) {
	 temp *= -1;
	 }
	 if (localDiff < temp) {
	 localDiff = temp;
	 }
	 }
	 }
      */
    for (int k = 0; k < numWorkers; k++){
      first = my_arg->j_init[k];
      last = my_arg->j_fin[k];

      ix = first*gridSize;
      for (i = first; i <= last; i++) {
	//int ix = i*gridSize;
	for (j = my_arg->j_init[myid]; j <= my_arg->j_fin[myid]; j++) {
	  grid1[ix+j] = (grid2[ix-gridSize+j] + grid2[ix+gridSize+j] + 
			 grid2[ix+j-1] + grid2[ix+j+1]) * 0.25;
	  temp = grid1[ix+j] - grid2[ix+j];
	  if (temp < 0) {
	    temp *= -1;
	  }
	  if (localDiff < temp) {
	    localDiff = temp;
	  }
	}
	ix+=gridSize;
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
  for (i = 1; i < gridSize-1; i++)
    for (j = 1; j < gridSize-1; j++) {
      grid1[IDX(gridSize,i,j)] = 0.0;
      grid2[IDX(gridSize,i,j)] = 0.0;
    }
  for (i = 0; i < gridSize; i++) {
    grid1[IDX(gridSize,i,0)] = 1.0;
    grid1[IDX(gridSize,i,gridSize-1)] = 1.0;
    grid2[IDX(gridSize,i,0)] = 1.0;
    grid2[IDX(gridSize,i,gridSize-1)] = 1.0;
  }
  for (j = 0; j < gridSize; j++) {
    grid1[IDX(gridSize,0,j)] = 1.0;
    grid2[IDX(gridSize,0,j)] = 1.0;
    grid1[IDX(gridSize,gridSize-1,j)] = 1.0;
    grid2[IDX(gridSize,gridSize-1,j)] = 1.0;
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
