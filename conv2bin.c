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

sample_t* read_sample(char *sample_fn, int *out_count)
{
	FILE *fp = fopen(sample_fn, "r");
	assert(fp!=NULL && "can't read sample file");
	sample_t *samples = malloc(sizeof(*samples) * 160000);
	assert(samples != NULL && "alloc error");
	int count=0;
	char *p = samples[0].name;
	size_t size = 64;
	while(getline(&p, &size, fp)>0)
	{
		strsep(&p, "\n");
		++count;
		p = samples[count].name;
	}
	samples = realloc(samples, sizeof(*samples) * count);
	assert(samples != NULL && "alloc error");
	qsort(samples, count, sizeof(*samples), strcmp_v);
	fclose(fp);

	*out_count = count;
	return samples;

}

void convert(char * in_fn, char * out_fn, sample_t* samples, int sample_count)
{
	char buff[1000];
	char * tok1, * tok2, *tok3, *next;

	item_t item={0};
	item.id1 = 0;
	

	gzFile in = gzopen(in_fn, "r");
	assert(in != Z_NULL && "gzopen in_fn failed");

	FILE * out = fopen(out_fn, "wb");
	assert(out != Z_NULL && "fopen out_fn failed");

	while(gzgets(in, buff, 1000)!=NULL)
	{
		next = buff;
		tok1 = strsep(&next, ":\t\n");
		assert(tok1 != NULL);
		tok2 = strsep(&next, ":\t\n");
		assert(tok2 != NULL);
		tok3 = strsep(&next, ":\t\n");
		assert(tok3 != NULL);
		while(strcmp_i(tok1, samples[item.id1].name)!= 0) {
			++item.id1;
			if(item.id1 >= sample_count) item.id1 = 0;
		}
		while(strcmp_i(tok2, samples[item.id2].name)!= 0) { 
			++item.id2;
			if(item.id2 >= sample_count) item.id2 = item.id1;
		}

		item.cM = strtof(tok3, NULL);

		/* write to output file */
		fwrite(&item, sizeof(item), 1, out);
	}


	gzclose(in);
	fclose(out);
}

int main(int argc, char **argv)
{
	char * usage = "Usage: ./conv2bin sample_list.txt in.gz out.gz";
	if (argc != 4)
	{
		fprintf(stderr, "%s\n", usage);
		exit(1);
	}
	int sample_count;
	sample_t * samples = read_sample(argv[1], &sample_count);
	convert(argv[2], argv[3], samples, sample_count);

	return 0;
}
