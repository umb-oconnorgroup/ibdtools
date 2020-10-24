#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

/* ----------------------------------------------------
 * vector of int
 * ---------------------------------------------------
 */
typedef struct {
    int *data;
    size_t size;
    size_t capacity;
} vint_t;

void
vint_alloc(vint_t *self, size_t size)
{
    assert(size >= 0);
    self->size = size;
    self->capacity = size;
    if (size > 0) {
        self->data = malloc(sizeof(*self->data) * size);
    } else {
        self->data = NULL;
    }
}

void
vint_free(vint_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vint_initialize(vint_t *self, int init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vint_add_element(vint_t *self, int value)
{
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = realloc(self->data, sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
}

/* ----------------------------------------------------
 * vector of double
 * ---------------------------------------------------
 */
typedef struct {
    double *data;
    size_t size;
    size_t capacity;
} vdouble_t;

void
vdouble_alloc(vdouble_t *self, size_t size)
{
    assert(size >= 0);
    self->size = size;
    self->capacity = size;
    if (size > 0) {
        self->data = malloc(sizeof(*self->data) * size);
    } else {
        self->data = NULL;
    }
}

void
vdouble_free(vdouble_t *self)
{
    if (self && self->data) {
        free(self->data);
        self->data = NULL;
    }
}

void
vdouble_init(vdouble_t *self, double init_value)
{
    for (int i = 0; i < self->size; i++) {
        self->data[i] = init_value;
    }
}

extern inline void
vdouble_add_element(vdouble_t *self, double value)
{
    if (self->size == self->capacity) {
        if (self->capacity == 0) {
            self->capacity = 1;
        } else {
            self->capacity *= 2;
        }
        self->data = realloc(self->data, sizeof(*self->data) * self->capacity);
        assert(self->data != NULL);
    }
    self->size += 1;
    self->data[self->size - 1] = value;
}

void
vdouble_fprintf(vdouble_t *self, FILE *fp)
{
    assert(fp != NULL);
    for (size_t i = 0; i < self->size; i++) {
        fprintf(fp, "%ld\t%g\n", i, self->data[i]);
    }
}

void
vdouble_fprintf_range(vdouble_t *self, FILE *fp, size_t start, size_t end)
{
    assert(fp != NULL);
    assert(start <= end);
    assert(end < self->size);
    for (size_t i = start; i <= end; i++) {
        fprintf(fp, "%ld\t%g\n", i, self->data[i]);
    }
}


