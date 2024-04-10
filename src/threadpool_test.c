#include <threadpool.h>

#include <pthread.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define NTHREAD 4

typedef struct job_arg_t {
  pthread_mutex_t mutex;
  unsigned int ncall;
} job_arg;

void perfect_job(void *arg) {
  job_arg *jarg = (job_arg *) arg;

  // Sleep 2 second
  sleep(2);

  // Increase the number of call performed
  pthread_mutex_lock(&jarg->mutex);
  jarg->ncall++;
  pthread_mutex_unlock(&jarg->mutex);
}

int check_state(threadpool_t *pool, job_arg *jarg, const int estatus, const unsigned int enthread, const unsigned int enactive, const unsigned int enjobs, const unsigned int encall) {
  const int status = threadpool_status(pool);
  unsigned int nthread, nactive, njobs;

  pthread_mutex_lock(&jarg->mutex);

  if(status != estatus) {
    fprintf(stderr, "Thread pool is in a wrong state (%04X instead of %04X)\n", status, estatus);
    threadpool_destroy(pool);
    pthread_mutex_unlock(&jarg->mutex);
    return -1;
  }

  if(threadpool_stats(pool, &nthread, &nactive, &njobs)) {
    fprintf(stderr, "Error while retrieving the pool statistics");
    pthread_mutex_unlock(&jarg->mutex);
    return -1;
  }

  if(nthread != enthread || nactive != enactive || njobs != enjobs) {
    fprintf(stderr, "The returned statistics are false (nthread=%u instead of %u, nactive=%u instead of %u, njobs=%u instead of %u)\n",
                    nthread, enthread, nactive, enactive, njobs, enjobs);
    pthread_mutex_unlock(&jarg->mutex);
    threadpool_destroy(pool);
    return -1;
  }

  if(jarg->ncall != encall) {
    fprintf(stderr, "The thread pool function was called %u times instead of %u\n", jarg->ncall, encall);
    pthread_mutex_unlock(&jarg->mutex);
    threadpool_destroy(pool);
    return -1;
  }

  pthread_mutex_unlock(&jarg->mutex);

  return 0;
}

int main(const int argc, const char **argv) {
  job_arg jarg = {PTHREAD_MUTEX_INITIALIZER, 0};
  threadpool_t *pool;

  // Create the thread pool
  pool = threadpool_create(NTHREAD);
  sleep(1); // Wait for the threads to be initialized

  // Check the status of the pool while running
  if(check_state(pool, &jarg, 0, NTHREAD, 0, 0, 0))
    return -1;

  // Add a job for each thread
  for(int i = 0; i < NTHREAD; i++)
    threadpool_addjob(pool, perfect_job, &jarg, NULL);
  sleep(1); // Wait for the threads to run

  // Check the status of the pool while running
  if(check_state(pool, &jarg, THREADPOOL_STATUS_RUNNING, NTHREAD, NTHREAD, NTHREAD, 0))
    return -1;

  // Wait for all jobs to finish
  threadpool_wait(pool);

  // Check the status of the pool after run
  if(check_state(pool, &jarg, 0, NTHREAD, 0, 0, NTHREAD))
    return -1;

  // Reset the number of call to 0
  jarg.ncall = 0;

  // Add 2*NTHREAD jobs
  for(int i = 0; i < NTHREAD; i++) {
    threadpool_addjob(pool, perfect_job, &jarg, NULL);
    threadpool_addjob(pool, perfect_job, &jarg, NULL);
  }
  sleep(1); // Wait for the threads to run

  // Cancel all threads, but still wait that the running jobs finished
  threadpool_cancel(pool, 0);

  // Check the status of the pool after cancelling
  if(check_state(pool, &jarg, 0, 0, 0, NTHREAD, NTHREAD))
    return -1;

  // Add NTHREAD threads in order to finish the job
  for(int i = 0; i < NTHREAD; i++)
    threadpool_addthread(pool);
  sleep(1); // Wait for the threads to run

  // Wait for all jobs to finish
  threadpool_wait(pool);

  // Check the status of the pool after taking back a cancelled job
  if(check_state(pool, &jarg, 0, NTHREAD, 0, 0, 2*NTHREAD))
    return -1;

  // Reset the number of call to 0
  jarg.ncall = 0;

  // Add NTHREAD+1 jobs
  for(int i = 0; i <= NTHREAD; i++)
    threadpool_addjob(pool, perfect_job, &jarg, NULL);
  sleep(1);

  // Cancel all threads, immediately
  threadpool_cancel(pool, -1);

  // Check the status of the pool after a 'brute' cancellation
  if(check_state(pool, &jarg, 0, 0, 0, 1, 0))
    return -1;

  // Add one thread
  threadpool_addthread(pool);
  sleep(1); // Sleep 1 seconds in order to be sure that the threads runs

  // Remove the thread, but wait for the whole jobs to be finished
  threadpool_removethread(pool, 1);

  // Wait for all jobs to finish
  threadpool_wait(pool);

  // Check the final status of the pool
  if(check_state(pool, &jarg, 0, 0, 0, 0, 1))
    return -1;

  // Destroy the pool
  threadpool_destroy(pool);

  sleep(2);
}
