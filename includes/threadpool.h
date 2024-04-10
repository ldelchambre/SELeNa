#ifndef _THREADPOOL_H_
#define _THREADPOOL_H_

/**
 * Thread pool structure and type definition
 */
typedef struct threadpool threadpool_t;

/**
 * Create a thread pool initially containing nthread threads
 *
 * @param[in] nthread The initial number of threads within the pool
 * @return A pointer to the allocated pool upon success, NULL upon failure.
 */
threadpool_t *threadpool_create(const unsigned int nthread);

/**
 * Add a job to a pool
 *
 * @param[in] pool The pool where to add the job
 * @param[in] fct The function to execute (if NULL, a job will be added that do nothing)
 * @param[in] arg User provided pointer to pass to the function to execute
 * @param[in] freearg A function designed the free arg once the job is done or cancelled (if NULL, no function will be called)
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_addjob(threadpool_t *pool, void (*fct)(void *), void *arg, void (*freearg)(void *));

/**
 * Add a new thread to the given pool.
 *
 * @param[in] pool The pool where to add a thread
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_addthread(threadpool_t *pool);

/**
 * Remove a thread from the pool.
 *
 * @param[in] pool The pool from which to remove the thread
 * @param[in] wait If non-zero, then wait for all the already-submitted jobs to finish before cancelling the thread.
 *
 *  Note: if the pool currently contains no threads, then the removal will be deferred to the moment
 *        the pool will contain at least one thread.
 *
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_removethread(threadpool_t *pool, const unsigned int wait);

/**
 * Bit-wise mask saying that the pool is currently running some task(s)
 */
#define THREADPOOL_STATUS_RUNNING    0x0001

/**
 * Bit-wise mask saying that the pool is currently cancelling all its threads
 */
#define THREADPOOL_STATUS_CANCELLING 0x0002

/**
 * Bit-wise mask saying that the pool is currently being destroyed
 */
#define THREADPOOL_STATUS_DESTROYING 0x0004

/**
 * Bit-wise mask saying that the pool is in an error state (i.e. the sole
 * permitted operation is threadpool_destroy)
 */
#define THREADPOOL_STATUS_ERROR      0x0008

/**
 * Get a bit-wise mask providing the current status of the pool.
 *
 * @param[in] pool The pool for which to retrieve the status
 * @return A bit-wise mask providing the current status of the pool
 */
int threadpool_status(const threadpool_t *pool);

/**
 * Get some basic statistics about the pool
 *
 * @param[in] pool The thread pool
 * @param[in] nthread If not NULL, return the current number of threads withi the pool
 * @param[in] nactive If not NULL, return the current number of active threads withi the pool
 * @param[in] njobs If not NULL, return the number of jobs within the pool (i.e. those being executed plus those being waiting)
 *
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_stats(threadpool_t *pool, unsigned int *nthread, unsigned int *nactive, unsigned int *njobs);

/**
 * Wait for all jobs to be finished.
 *
 * @param[in] pool The pool to wait for
 *
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_wait(threadpool_t *pool);

/**
 * Cancel all threads from a pool.
 *
 * @param[in] pool The pool
 * @param[in] wait The cancellation schedule for the threads
 *
 * Cancellation schedule is as follow:
 *   - If wait < 0: Cancel all threads immediately, without waiting for the running jobs to be finised
 *   - If wait = 0: Cancel all threads after the execution of the running jobs
 *   - If wait > 0: Cancel all threads after all currently submitted jobs are finished
 *
 * @return 0 upon success; -1 is the pool is in an error state or an error code upon other failure
 */
int threadpool_cancel(threadpool_t *pool, const int wait);

/**
 * Destroy a pool.
 *
 * @param[in] pool The pool to destroy
 */
void threadpool_destroy(threadpool_t *pool);

#endif
