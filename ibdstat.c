#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#define MEM_UNIT 10000

struct IBD
{
	double *arr;
	size_t size;
	size_t max;
};

void
IBD_init(struct IBD *ibd)
{
	ibd->max = MEM_UNIT;
	ibd->arr = (double *)malloc(sizeof(double) * ibd->max);
	ibd->size = 0;
	assert(ibd->arr != NULL);
}

void
IBD_add(struct IBD* ibd, double cM)
{
	if(ibd->size + 10 > ibd->max)
	{
		ibd->max += MEM_UNIT;
		ibd->arr = realloc(ibd->arr, ibd->max * sizeof(double));
		assert(ibd->arr != NULL);
	}
	ibd->arr[ibd->size] = cM;
	ibd->size += 1;

//	printf("size= %ld, max= %ld\n", ibd->size, ibd->max);
}

int 
IBD_compare(const void * d1, const void * d2)
{
	return *((double *)d1) > *((double *)d2) ? 1 : -1;
}

void 
IBD_sort(struct IBD* ibd)
{
	qsort((void*)ibd->arr, ibd->size, sizeof(double), IBD_compare);
}


// should not be called before sorting
void
IBD_print_quantiles(struct IBD *ibd)
{
	int quantile_size = 1000;

	size_t index=0;
	for(size_t i=0; i< quantile_size; i++)
	{
		index = (size_t)( 1.0 * i / quantile_size * ibd->size) ;
		if(index == ibd->size) index = ibd->size - 1;
		fprintf(stderr, "%g\t%g\n", 1.0 * i / quantile_size, ibd->arr[index]);
	}
}


// should not be called before sorting
void
IBD_print_hist(struct IBD *ibd)
{
	size_t index=0, cnt_in_range=0;
	double start=0, step=1, end=step;
	while(index < ibd->size)
	{
		if(ibd->arr[index] < end) {
			cnt_in_range++;
			index ++;
		}
		else{
			// print current range and its counts.
			fprintf(stdout, "%g\t%g\t%ld\t%g\n", start, end, cnt_in_range, 1.0 * cnt_in_range / ibd->size);
			// prepare for next range
			start = end;
			end = end + step;
			cnt_in_range = 0;
			// index no change put it back and check again.
		}
	}
	fprintf(stdout, "%g\t%g\t%ld\t%g\n", start, end, cnt_in_range, 1.0 * cnt_in_range / ibd->size);
}

void 
IBD_print(struct IBD *ibd)
{
	printf("# -------------- hist\n");
	IBD_print_hist(ibd);
	printf("# -------------- quantiles\n");
	IBD_print_quantiles(ibd);
}

void IBD_free(struct IBD* ibd)
{
	if(ibd->arr) free(ibd->arr);
	ibd->arr = NULL;
}


int main()
{
	char buff[1000];
	char *p = buff, *pChar = NULL;
	size_t size = 999;

	double cM=0, cM_min=100, cM_max=0, cM_sum=0, cM_avg=0;
	size_t cM_count=0;

	int i=0, colCount=0;
	struct IBD ibd;

	IBD_init(&ibd);
	while( getline(&p, &size, stdin) > 0)
	{
		cM_count++;	
		// find the right column
		for(pChar=p, colCount=1; *pChar != '\0'; pChar++)
		{
			if (*pChar == ' ' || *pChar == '\t') ++colCount;
			if (colCount == 8) break;
		}
		assert(colCount==8);
		++pChar;
		cM = strtod(pChar, NULL);

		cM_sum += cM;
		cM_min = cM_min > cM ? cM : cM_min;
		cM_max = cM_max > cM ? cM_max : cM;

		IBD_add(&ibd, cM);
	}

	cM_avg = cM_sum / cM_count;

	IBD_sort(&ibd);

	printf("# n=%ld, avg=%lf, min=%lf, max=%lf\n", cM_count, cM_avg, cM_min, cM_max);

	IBD_print(&ibd);

	IBD_free(&ibd);
	return 0;
}
