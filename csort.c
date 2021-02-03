/*
 * Key points:
 * 	1. take advantage of the sorted id order. Use linear search rather than others
 * 	2. Use inline version of strcmp_i instead of library strcmp function.
 */
#define _GNU_SOURCE
#include<search.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<zlib.h>
#include<string.h>
#include<stdint.h>
#define TOTAL_CM 3545.83

typedef struct {
	uint32_t id1, id2;
	float cM; 
} item_t;

typedef struct {
	char name[96];
} sample_t;

/* inline version of strcmp */
static inline int
strcmp_i(const char *a, const char *b)
{
    while (1) {
        if (*a != *b)
            return *a - *b;
        else {
            if (*a == '\0')
                return 0;
            ++a;
            ++b;
        }
    }
}

int strcmp_v(const void *a, const void *b)
{
	return strcmp_i(a, b);
}

void read_sample()
{
	sample_t *samples = malloc(sizeof(*samples) * 160000);
	assert(samples != NULL && "alloc error");
	int count=0;
	char *p = samples[0].name;
	size_t size = 64;
	while(getline(&p, &size, stdin)>0)
	{
		strsep(&p, "\n");
		++count;
		p = samples[count].name;
	}
	samples = realloc(samples, sizeof(*samples) * count);
	assert(samples != NULL && "alloc error");
	qsort(samples, count, sizeof(*samples), strcmp_v);
	for(int i=0; i<count; i++){
		printf("%s\n", samples[i].name);
	}
}


int main(int argc, char **argv)
{
	read_sample();

	return 0;
}

