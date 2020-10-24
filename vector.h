#include<stdlib.h>
#include<stdio.h>

typedef struct {
    int *data;
    size_t size;
    size_t capacity;
} vint_t;

void vint_alloc(vint_t *self, size_t size);
void vint_free(vint_t *self);
void vint_initialize(vint_t *self, int init_value);
void vint_add_element(vint_t *self, int value);

typedef struct {
    double *data;
    size_t size;
    size_t capacity;
} vdouble_t;

void vdouble_alloc(vdouble_t *self, size_t size);
void vdouble_free(vdouble_t *self);
void vdouble_init(vdouble_t *self, double init_value);
void vdouble_add_element(vdouble_t *self, double value);
void vdouble_fprintf(vdouble_t *self, FILE *fp);
void vdouble_fprintf_range(vdouble_t *self, FILE *fp, size_t start, size_t end);
