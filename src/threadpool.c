#include <threadpool.h>

#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Comment in order to disable the debug mode
 */
// #define _THREADPOOL_DEBUG_MODE

/**
 * In debug mode, _THREADPOOL_DEBUG will print a debug message
 */
#ifdef _THREADPOOL_DEBUG_MODE
#define _THREADPOOL_DEBUG(format,...) do{fprintf(stdout, format, ##__VA_ARGS__);fflush(stdout);}while(0)
#else
#define _THREADPOOL_DEBUG(format,...)
#endif

/**
 * Common error handling within the thread pool
 */
#define _THREADPOOL_ERROR(err,action, ...)\
do {\
  int errsv = err; \
  fprintf(stderr, "Threadpool - Error while "action": ", ##__VA_ARGS__); \
  fprintf(stderr, "%s\n", strerror(errsv)); \
  fflush(stderr);\
} while(0)

/**
 * Job node
 */
typedef struct thread_job_t {
  void (*fct)(void *); // The function to execute (if NULL nothing is done)
  void *arg; // The argument to the function
  void (*freearg)(void *); // The function to call in order to free arg (NULL if none)
  struct thread_job_t *prev; // The previous jobs
  struct thread_job_t *next; // The next job
} thread_job;

/**
 * Thread node
 */
typedef struct thread_queue_t {
  pthread_t id; // The thread identifier
  thread_job *job; // The job currently executed by this thread (NULL if thread is currently waiting for a job)
  threadpool_t *pool; // The pool this thread belongs to
  struct thread_queue_t *prev; // The previous thread within the pool
  struct thread_queue_t *next; // The next thread within the pool
} thread_queue;

/**
 * Double-ended queue of jobs within a pool
 */
typedef struct thread_job_queue_t {
  thread_job *front; // The front end of the queue
  thread_job *back; // The back end of the queue
} thread_job_queue;

/**
 * Thread pool structure
 */
struct threadpool {
  pthread_mutex_t lock; // Lock on the pool (allow thread synchronization)

  thread_queue *threads; // The list of threads
  thread_job_queue jobs_waiting; // The jobs that are currently waiting to be executed
  thread_job_queue jobs_running; // The jobs being currently executed

  pthread_cond_t job_available; // Condition signalling that a new job has been added to the pool
  pthread_cond_t jobs_finished; // Condition signalling that all jobs were finished
  pthread_cond_t threads_finished; // Condition signalling the end of all threads

  int status; // The current status of the pool
};

/**
 * Add a job node into the double ended queue of jobs specified by queue.
 * If wait is 0, the job node is inserted in front of the queue, otherwise
 * at the end of the queue.
 * Note: This function is THREAD-UNSAFE.
 */
void _threadpool_addjob_queue(thread_job_queue *queue, thread_job *job, const int wait) {
  if(queue->front == NULL) { // The list is currently empty
    queue->front = queue->back = job;
    job->next = NULL;
    job->prev = NULL;
  } else if(wait) { // Add this job at the end of the list, if we are allowed to wait
    queue->back->next = job;
    job->prev = queue->back;
    job->next = NULL;
    queue->back = job;
  } else { // Add this job in front of the list if we cannot wait
    queue->front->prev = job;
    job->prev = NULL;
    job->next = queue->front;
    queue->front = job;
  }
}

/**
 * Remove the job node from the specified queue.
 * Note that the job MUST be contained within the queue.
 * Note: This function is THREAD-UNSAFE.
 */
void _threadpool_removejob_queue(thread_job_queue *queue, const thread_job *job) {
  if(job->prev)
    job->prev->next = job->next;
  else
    queue->front = job->next;
  if(job->next)
    job->next->prev = job->prev;
  else
    queue->back = job->prev;
}

/**
 * Add the calling thread into the list of threads of the given pool.
 * @return a pointer to the thread node upon success or NULL upon failure
 * Note: This function is THREAD-UNSAFE.
 */
thread_queue *_threadpool_addthread(threadpool_t *pool) {
  thread_queue *thread;

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Adding this thread into the pool\n", pool, (unsigned) pthread_self(), pool->status);

  // Allocate the thread
  thread = (thread_queue *) malloc(sizeof(thread_queue));
  if(thread == NULL) {
    _THREADPOOL_ERROR(ENOMEM, "allocating a new thread");
    return NULL;
  }

  // Initialize the thread
  thread->id = pthread_self();
  thread->job = NULL;
  thread->pool = pool;

  // Add this thread into the pool
  thread->prev = NULL;
  thread->next = pool->threads;
  pool->threads = thread;
  if(thread->next)
    thread->next->prev = thread;

  return thread;
}

/**
 * Remove the given thread from the list of threads associated with the
 * pool it belongs to.
 * Note: This function is THREAD-UNSAFE.
 */
void _threadpool_removethread(thread_queue *thread) {
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Removing thread %012X from the pool\n", thread->pool, (unsigned) pthread_self(), thread->pool->status, (unsigned)thread->id);

  if(thread->prev)
    thread->prev->next = thread->next;
  else
    thread->pool->threads = thread->next;
  if(thread->next)
    thread->next->prev = thread->prev;
  free(thread);
}

/**
 * Add a job into a pool.
 * If wait is 0, then the job will be the next one to be executed, otherwise
 * this will be the last one.
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 * Note: This function is THREAD-UNSAFE.
 */
int _threadpool_addjob(threadpool_t *pool, void (*fct)(void *), void *arg, void (*freearg)(void *), const int wait) {
  thread_job *job;
  int err;

  // If we are in an error state, no new job can be added
  if(pool->status & THREADPOOL_STATUS_ERROR) {
    if(freearg) freearg(arg);
    return -1;
  }

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Adding a job (fct=%p, arg=%p, freearg=%p)\n", pool, (unsigned) pthread_self(), pool->status, fct, arg, freearg);

  // Allocate the job
  job = (thread_job *) malloc(sizeof(thread_job)); // XXX If many short jobs are added to the queue, then malloc constitute a bottleneck
  if(job == NULL) {
    _THREADPOOL_ERROR(ENOMEM, "allocating a new job");
    if(freearg) freearg(arg);
    return ENOMEM;
  }

  // Initialize the job
  job->fct = fct;
  job->arg = arg;
  job->freearg = freearg;

  // Add this jobs to the job waiting list
  _threadpool_addjob_queue(&pool->jobs_waiting, job, wait);

  // Signal the pool that a new job is available
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Signallin a new job (fct=%p, arg=%p, freearg=%p)\n", pool, (unsigned) pthread_self(), pool->status, fct, arg, freearg);
  err = pthread_cond_signal(&pool->job_available);
  if(err) {
    _THREADPOOL_ERROR(err, "signaling a new job");
    _threadpool_removejob_queue(&pool->jobs_waiting, job);
    if(freearg) freearg(arg); 
    free(job);
    return err;
  }

  return 0;
}

/**
 * Cancel all thread from a pool.
 * If wait is:
 *   - wait < 0 : We are pressed, so cancel all threads immediately
 *   - wait == 0: Wait for the jobs to end their execution before cancelling
 *   - wait > 0 : Wait for the whole jobs to be executed before cancelling
 * @return 0 upon success or an error code upon failure
 * Note: This function is THREAD-UNSAFE.
 */
int _threadpool_cancel(threadpool_t *pool, const int wait) {
  thread_queue *thread;
  int err;

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Cancelling all threads (wait=%d)\n", pool, (unsigned) pthread_self(), pool->status, wait);

  // Set this pool in a cancellation state
  pool->status |= THREADPOOL_STATUS_CANCELLING;

  // Cancel all threads
  thread = pool->threads;
  while(thread) {

    if(wait < 0) // If wait < 0, 'immediately' destroy the threads
      err = pthread_cancel(thread->id);
    else // Otherwise, add exit jobs in order to cancel all threads
      err = _threadpool_addjob(pool, pthread_exit, NULL, NULL, (wait > 0));

    // Upon error
    if(err) {
      _THREADPOOL_ERROR(err, "cancelling thread");
      // If we are not in a destroy state, then exit on error
      // otherwise keep cancelling threads
      if((pool->status & THREADPOOL_STATUS_DESTROYING) == 0) {
        pool->status &= ~THREADPOOL_STATUS_CANCELLING;
        return err;
      }
    }

    // Go to next thread
    thread = thread->next;
  }

  // Wait for the whole threads to ends
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Waiting for all threads to finish\n", pool, (unsigned) pthread_self(), pool->status);
  while(pool->threads) {
    if(wait < 0) { // If wait == -1, use a timed condition wait
      struct timespec t;
      clock_gettime(CLOCK_REALTIME, &t);
      t.tv_sec += 60; // Set a timeout of 1 minute
      err = pthread_cond_timedwait(&pool->threads_finished, &pool->lock, &t);
    } else // Otherwise, use a classical wait
      err = pthread_cond_wait(&pool->threads_finished, &pool->lock);
    if(err)
      _THREADPOOL_ERROR(err, "waiting for all threads to finish");
  }
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - All threads finished\n", pool, (unsigned) pthread_self(), pool->status);

  // We are no more in a cancellation state
  pool->status &= ~THREADPOOL_STATUS_CANCELLING;

  return err;
}

/**
 * Destroy a pool (i.e. cancel threads and flush the jobs lists)
 */
void _threadpool_destroy(threadpool_t *pool) {
  thread_queue *thread;
  thread_job *job;

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Destroying the pool\n", pool, (unsigned) pthread_self(), pool->status);

  // Set this pool in a destroy state
  pool->status |= THREADPOOL_STATUS_DESTROYING;

  // Cancel all threads immediately
  _threadpool_cancel(pool, -1);

  // Flush the waiting list
  while(pool->jobs_waiting.front) {
    job = pool->jobs_waiting.front;
    pool->jobs_waiting.front = job->next;
    if(job->freearg) job->freearg(job->arg); 
    free(job);
  }

  // Flush the threads list (this should be empty but just in case of)
  while(pool->threads) {
    thread = pool->threads;
    pool->threads = thread->next;
    free(thread);
  }

  // Flush the running list (this should be empty but just in case of)
  while(pool->jobs_running.front) {
    job = pool->jobs_running.front;
    pool->jobs_running.front = job->next;
    if(job->freearg) job->freearg(job->arg); 
    free(job);
  }

  // Notify all waiters that all jobs were finished
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Notifying the end of all jobs\n", pool, (unsigned) pthread_self(), pool->status);
  pthread_cond_broadcast(&pool->jobs_finished);
}

/**
 * Set this pool in an error state, then exit.
 */
void _threadpool_error(threadpool_t *pool) {
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Switching to error state", pool, (unsigned) pthread_self(), pool->status);

  // Set this pool in an error state (i.e. don't accept, nor run any further jobs)
  pool->status |= THREADPOOL_STATUS_ERROR;

  // Exit this thread
  pthread_exit(NULL);
}

/**
 * Lock the pool
 * @return 0 upon success or an error code upon failure
 */
int _threadpool_lock(threadpool_t *pool) {
  const int err = pthread_mutex_lock(&pool->lock);
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Pool locked\n", pool, (unsigned) pthread_self(), pool->status);
  if(err)
    _THREADPOOL_ERROR(err, "acquiring the pool lock");
  return err;
}

/**
 * Unlock the pool
 * @return 0 upon success or an error code upon failure
 */
int _threadpool_unlock(threadpool_t *pool) {
  const int err = pthread_mutex_unlock(&pool->lock);
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Pool unlocked\n", pool, (unsigned) pthread_self(), pool->status);
  if(err)
    _THREADPOOL_ERROR(err, "releasing the pool lock");
  return err;
}

/**
 * Clean the ressources handled by a thread.
 * Note: The pool lock MUST be acquired by this thread before calling
 *       _threadpool_threadcleanup. The lock will then be released.
 */
void _threadpool_threadcleanup(void *arg) {
  thread_queue *thread_node  = (thread_queue *) arg;
  threadpool_t *pool = thread_node->pool;
  int err;

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Cleaning up the thread ressources\n", pool, (unsigned) pthread_self(), pool->status);

  // Remove this thread from the pool
  _threadpool_removethread(thread_node);

  // If no more thread exist, signal all waiters
  if (pool->threads == NULL) {
    _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Notifying the end of all threads\n", pool, (unsigned) pthread_self(), pool->status);
    err = pthread_cond_signal(&pool->threads_finished);
    if(err) {
      _THREADPOOL_ERROR(err, "signaling the end of all threads");
      _threadpool_error(pool);
    }
  }

  // Release the pool lock
  if(_threadpool_unlock(pool))
    _threadpool_error(pool);
}

/**
 * Clean the ressources handled by a job.
 */
void _threadpool_jobcleanup(void *arg) {
  thread_queue *thread_node  = (thread_queue *) arg;
  threadpool_t *pool = thread_node->pool;
  int err;

  // Acquire the pool lock
  if(_threadpool_lock(pool))
    _threadpool_error(pool);

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Cleaning up the job ressources (fct=%p, arg=%p)\n", pool, (unsigned) pthread_self(), pool->status, thread_node->job->fct, thread_node->job->arg);

  // Remove this job from the running list
  _threadpool_removejob_queue(&pool->jobs_running, thread_node->job);

  // Delete this jobs
  if(thread_node->job->freearg)
    thread_node->job->freearg(thread_node->job->arg);
  free(thread_node->job);

  // The thread that were running this job don't have a job anymore
  thread_node->job = NULL;

  // If all running jobs are done
  if(pool->jobs_running.front == NULL) {
    // Thread pool is in a wait state now
    pool->status &= ~THREADPOOL_STATUS_RUNNING;

    // If all jobs are done
    if(pool->jobs_waiting.front == NULL) {
      // Notify the waiters
      _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Notifying the end of all jobs\n", pool, (unsigned) pthread_self(), pool->status);
      err = pthread_cond_broadcast(&pool->jobs_finished);
      if(err) {
        _THREADPOOL_ERROR(err, "broadcasting the end of all jobs");
        _threadpool_error(pool);
        return;
      }
    }
  }
}

/**
 * Thread Pool main loop
 */
void *_threadpool_main(void *arg) {
  threadpool_t *pool = (threadpool_t *) arg;
  thread_queue *thread_node;
  sigset_t pset; // The set of all signals to catch (inherited from the launching thread)
  int err;

  // Defer any cancellation of this thread such that we can take the time
  // to correctly set up the thread within the pool
  pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &err);

  // Get the signal mask of the calling thread
  pthread_sigmask(0, NULL, &pset);

  // Acquire the pool's lock
  if(_threadpool_lock(pool))
    _threadpool_error(pool);

  // Add me to the pool
  thread_node = _threadpool_addthread(pool);
  if(thread_node == NULL)
    _threadpool_error(pool);

  // Push the thread cleanup procedure
  pthread_cleanup_push(_threadpool_threadcleanup, thread_node);

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Running the main loop\n", pool, (unsigned) pthread_self(), pool->status);

  // thread pool main loop
  while((pool->status & THREADPOOL_STATUS_ERROR) == 0) {
    // Reset signal mask and re-enable deferred cancellation requests
    // Note: this stands here because these might have been modified by
    //       the job function
    pthread_sigmask(SIG_SETMASK, &pset, NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, &err);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &err);

    // Wait for some jobs to become available
    _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Waiting for jobs to become available\n", pool, (unsigned) pthread_self(), pool->status);
    while(pool->jobs_waiting.front == NULL) {
      err = pthread_cond_wait(&pool->job_available, &pool->lock);
      if(err) {
        _THREADPOOL_ERROR(err,"waiting for new jobs");
        _threadpool_error(pool);
      }
    }

    // If we are in an error state, then exit
    if(pool->status & THREADPOOL_STATUS_ERROR)
      pthread_exit(NULL);

    // Defer cancellation requests
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &err);

    // Assign the next job to this thread
    thread_node->job = pool->jobs_waiting.front;

    // Remove the job from the waiting list
    _threadpool_removejob_queue(&pool->jobs_waiting, thread_node->job);

    // Add it to the running list
    _threadpool_addjob_queue(&pool->jobs_running, thread_node->job, 0);

    // Thread pool is in a running state
    pool->status |= THREADPOOL_STATUS_RUNNING;

    // Release the pool's lock
    if(_threadpool_unlock(pool))
      _threadpool_error(pool);

    // Push the job cleanup procedure
    pthread_cleanup_push(_threadpool_jobcleanup, thread_node);

    // Re-enable thread cancellation
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &err);

    // Run the job
    _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Running job (fct=%p, arg=%p)\n", pool, (unsigned) pthread_self(), pool->status, thread_node->job->fct, thread_node->job->arg);
    if(thread_node->job->fct)
      thread_node->job->fct(thread_node->job->arg);

    // Remove job from the pool
    pthread_cleanup_pop(1);
  }

  // Remove and execute the thread cleanup procedure
  pthread_cleanup_pop(1);
  pthread_exit(NULL);
}

threadpool_t *threadpool_create(const unsigned int nthread) {
  threadpool_t *pool;
  int err;

  // Allocate the pool
  pool = (threadpool_t *) malloc(sizeof(threadpool_t));
  if(pool == NULL) {
    _THREADPOOL_ERROR(ENOMEM, "allocating the pool");
    return NULL;
  }

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x00) - Creating a pool with %u threads\n", pool, (unsigned) pthread_self(), nthread);

  // Initialize the thread pool
  pool->threads = NULL;
  pool->jobs_waiting.front = NULL;
  pool->jobs_waiting.back = NULL;
  pool->jobs_running.front = NULL;
  pool->jobs_running.back = NULL;
  pool->status = 0;
  if((err = pthread_mutex_init(&pool->lock, NULL))) {
    _THREADPOOL_ERROR(err, "initializing the pool mutex");
    free(pool);
    return NULL;
  }
  if((err = pthread_cond_init(&pool->job_available, NULL))) {
    _THREADPOOL_ERROR(err, "initializing the pool job synchronization");
    pthread_mutex_destroy(&pool->lock);free(pool);
    return NULL;
  }
  if((err = pthread_cond_init(&pool->jobs_finished, NULL))) {
    _THREADPOOL_ERROR(err, "initializing the pool wait synchronization");
    pthread_cond_destroy(&pool->job_available);pthread_mutex_destroy(&pool->lock);free(pool);
    return NULL;
  }
  if((err = pthread_cond_init(&pool->threads_finished, NULL))) {
    _THREADPOOL_ERROR(err, "initializing the pool threads synchronization");
    pthread_cond_destroy(&pool->jobs_finished);pthread_cond_destroy(&pool->job_available);pthread_mutex_destroy(&pool->lock);free(pool);
    return NULL;
  }

  // Launch all threads
  for(unsigned int i = 0; i < nthread; i++) {
    if(threadpool_addthread(pool)) {
      threadpool_destroy(pool);
      return NULL;
    }
  }

  return pool;
}

int threadpool_addthread(threadpool_t *pool) {
  pthread_attr_t attr;
  pthread_t id;
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR)
    return -1;

  // Initialize the thread attributes
  err = pthread_attr_init(&attr);
  if(err) {
    _THREADPOOL_ERROR(err, "initializing the thread attributes");
    return err;
  }

  // Created threads are not intended to be joinable
  err = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  if(err) {
    _THREADPOOL_ERROR(err, "setting the thread detach state");
    pthread_attr_destroy(&attr);
    return err;
  }

  // Launch the thread along with the provided attributes
  err = pthread_create(&id, &attr, _threadpool_main, pool);
  if(err) {
    _THREADPOOL_ERROR(err, "launching thread");
    pthread_attr_destroy(&attr);
    return err;
  }

  // Destroy the thread attributes
  err = pthread_attr_destroy(&attr);
  if(err) {
    _THREADPOOL_ERROR(err, "destroying the thread attributes");
    pthread_cancel(id);
    return err;
  }

  return 0;
}

int threadpool_removethread(threadpool_t *pool, const unsigned int wait) {
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR)
    return -1;

  // Acquire the pool's lock
  err = _threadpool_lock(pool);
  if(err)
    return err;

  // Add a job that will make a thread to exit
  err = _threadpool_addjob(pool, pthread_exit, NULL, NULL, wait);
  if(err) {
    _threadpool_unlock(pool);
    return err;
  }

  // Release the pool's lock
  err = _threadpool_unlock(pool);
  if(err)
    return err;

  return 0;
}

int threadpool_status(const threadpool_t *pool) {
  return pool->status;
}

int threadpool_stats(threadpool_t *pool, unsigned int *nthread, unsigned int *nactive, unsigned int *njobs) {
  unsigned int _nthread, _nactive, _njobs;
  thread_queue *thread;
  thread_job *job;
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR)
    return -1;

  // If we don't have to compute anything, quit
  if(nthread == NULL && nactive == NULL && njobs == NULL)
    return 0;

  // Acquire the pool's lock
  err = _threadpool_lock(pool);
  if(err)
    return err;

  // Count the number of threads
  _nthread = _nactive = 0;
  thread = pool->threads;
  while(thread) {
    _nthread++;
    _nactive += (thread->job != NULL);
    thread = thread->next;
  }

  // Count the number of waiting jobs
  _njobs = 0;
  job = pool->jobs_waiting.front;
  while(job) {
    _njobs++;
    job = job->next;
  }

  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Getting the pool's statistics (nthread=%u, nactive=%u, njobs=%u)\n", pool, (unsigned) pthread_self(), pool->status, _nthread, _nactive, _nactive+_njobs);

  // Release the pool's lock
  err = _threadpool_unlock(pool);
  if(err)
    return err;

  // Set the stats
  if(nthread) *nthread = _nthread;
  if(nactive) *nactive = _nactive;
  if(njobs) *njobs = _njobs + _nactive;

  return err;
}

int threadpool_addjob(threadpool_t *pool, void (*fct)(void *), void *arg, void (*freearg)(void *)) {
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR) {
    if(freearg) freearg(arg);
    return -1;
  }

  // Acquire the pool's lock
  err = _threadpool_lock(pool);
  if(err) {
    if(freearg) freearg(arg);
    return err;
  }

  // Create a new job
  err = _threadpool_addjob(pool, fct, arg, freearg, 1);
  if(err) {
    _threadpool_unlock(pool);
    return err;
  }

  // Release the pool's lock
  err = _threadpool_unlock(pool);
  if(err)
    return err;

  return 0;
}

int threadpool_wait(threadpool_t *pool) {
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR)
    return -1;

  // Acquire the pool's lock
  err = _threadpool_lock(pool);
  if(err)
    return err;

  // Wait for all jobs to finish
  while(pool->jobs_waiting.front != NULL || pool->jobs_running.front != NULL) {
    _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - Waiting for all jobs to finish\n", pool, (unsigned) pthread_self(), pool->status);
    err = pthread_cond_wait(&pool->jobs_finished, &pool->lock);
    if(err) {
      _THREADPOOL_ERROR(err,"waiting for all jobs to finish");
      _threadpool_unlock(pool);
      return err;
    }
  }
  _THREADPOOL_DEBUG("Threadpool %p (thread=%012X,status=0x%02X) - All jobs finished\n", pool, (unsigned) pthread_self(), pool->status);

  // Release the pool's lock
  err = _threadpool_unlock(pool);
  if(err)
    return err;

  return 0;
}

int threadpool_cancel(threadpool_t *pool, const int wait) {
  int err;

  // If we are in an error state, don't do anything
  if(pool->status & THREADPOOL_STATUS_ERROR)
    return -1;

  // Acquire the pool's lock
  err = _threadpool_lock(pool);
  if(err)
    return err;

  // Cancel all threads
  err = _threadpool_cancel(pool, wait);
  if(err) {
    _threadpool_unlock(pool);
    return err;
  }

  // Release the pool's lock
  err = _threadpool_unlock(pool);
  if(err)
    return err;

  return 0;
}

void threadpool_destroy(threadpool_t *pool) {
  if(pool == NULL) return;
  _threadpool_lock(pool); // Acquire the lock
  _threadpool_destroy(pool); // Destroy the pool

  // Destroy the synchronization variables
  pthread_cond_destroy(&pool->job_available);
  pthread_cond_destroy(&pool->jobs_finished);
  pthread_cond_destroy(&pool->threads_finished);
  _threadpool_unlock(pool); // Release the lock
  pthread_mutex_destroy(&pool->lock);

  // Delete this pool
  free(pool);
}
