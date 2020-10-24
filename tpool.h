#ifndef _tpool_h_
#define _tpool_h_

#include<pthread.h>

typedef struct {
    int num_threads;
    pthread_mutex_t mutex;
    pthread_t *threads;
    long running_id;
    long last_id;
    long step;
    void *(*job_func)(void *, long);
    void *parent_void;
}tpool_t;

void tpool_alloc(tpool_t *self, void* parent_void,  int num_threads);
void tpool_run(tpool_t *self, long thread_init_job_id, long thread_last_job_id, 
	long step, void *(*job_func)(void *, long));
void tpool_free(tpool_t *self);

#endif
