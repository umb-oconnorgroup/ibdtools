#include<pthread.h>
#include<stdlib.h>
#include<assert.h>

/* ----------------------------------------------------
 * struct for multithreading
 * ---------------------------------------------------
 */

typedef struct {
    int num_threads;
    pthread_mutex_t mutex;
    pthread_t *threads;
    long running_id;
    long last_id;
    long step;
    void *(*job_func)(void *, long);
    void *parent_void;
} tpool_t;

void
tpool_alloc(tpool_t *self, void* parent_void,  int num_threads)
{
    if (num_threads < 1)
        num_threads = 1;
    self->num_threads = num_threads;
    self->threads = malloc(sizeof(*self->threads) * num_threads);
    self->parent_void = parent_void;
    assert(self->threads != NULL);
}

// tpool_run call tpool_thread_func, and tpool_thread_func call parent job_func.
static inline void *
tpool_thread_func (void * self_void) 
{
    tpool_t *self = (tpool_t *) self_void;
    long start_job_id, cur_job_id;

    while (1) {
        // get starting job id
        pthread_mutex_lock(&self->mutex);
        start_job_id = self->running_id;
        if (start_job_id < self->last_id) {
            pthread_mutex_unlock(&self->mutex);
            return NULL;
        }
        self->running_id -= self->step;
        pthread_mutex_unlock(&self->mutex);

        // do the jobs: each thread do `step` jobs per batch
        for (long i = 0; i < self->step; i++) { 
		cur_job_id = start_job_id - i ;
		if (cur_job_id < self->last_id) {
			break;
		}
	 	self->job_func(self->parent_void, cur_job_id);
        }
    }
    return NULL;
}

void
tpool_run(tpool_t *self, long thread_init_job_id, long thread_last_job_id, 
	long step, void *(*job_func)(void *, long))
{
    self->job_func = job_func;
    self->running_id = thread_init_job_id;
    self->last_id = thread_last_job_id;
    self->step = step;

    // fprintf(stderr, "called: tpool_run\n");
    // fprintf(stderr, "num_thread: %d\n", self->num_threads);
    for (int i = 0; i < self->num_threads; i++) {
        pthread_create(self->threads + i, NULL, tpool_thread_func, (void *) self);
    }
    for (int i = 0; i < self->num_threads; i++) {
        pthread_join(self->threads[i], NULL);
    }
}


void
tpool_free(tpool_t *self)
{
    if (self->threads != NULL) {
        free(self->threads);
        self->threads = NULL;
        self->num_threads = 0;
    }
}
