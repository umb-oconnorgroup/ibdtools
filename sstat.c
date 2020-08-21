#include<stdio.h>
#include<string.h>
#include<assert.h>
#include<stdlib.h>

struct vector
{
	double *arr;
	size_t nmemb;
	size_t limit;
	size_t inc;
};

int cmp_double(const void * d1, const void * d2)
{
	double *x1 = (double *) d1;
	double *x2 = (double *) d2;
	if(*x1 == *x2) return 0;
	else if(*x1 > *x2) return 1;
	else return -1;
}

int main(int argc, char * argv[])
{
	char *linePtr=NULL;
	size_t lineSize=0;
	int nquan = 1000;
	int col=1;
	double bin_width = 1.0;
	double first_lower_bound=0;

	if(argc == 5) 
	{	col = atoi(argv[1]);
		nquan = atoi(argv[2]);
		first_lower_bound = strtod(argv[3], NULL);
		bin_width = strtod(argv[4], NULL);
	}
	else
	{
		fprintf(stderr, "\n\nUsage: cat tab_delim_file | \\\n"
				"           ./sstat <int_ith_col> <int_quantile> <dbl_first_lower_bound> <dbl_bin_width> 1>quantile.txt 2>histogram.txt \n\n");
		exit(1);
	}

	struct vector v;
	v.limit = 1000000;
	v.inc = 1000000;
	v.arr = (double *)malloc(sizeof(*v.arr) * v.limit);
	assert(v.arr!=NULL);
	v.nmemb = 0;

	while(getline(&linePtr, &lineSize, stdin)>0 && *linePtr != '\n')
	{
		char * token = strtok(linePtr, "\t ");
		for(int i = 2; i<= col; i++) token = strtok(NULL, "\t ");
		double value = strtod(token, NULL);
		
		if(v.nmemb + 10 > v.limit) 
		{
			v.limit += v.inc;
			v.arr = realloc(v.arr, sizeof(*v.arr) * v.limit);
			assert(v.arr!=NULL);
		}
		v.arr[v.nmemb] = value;
		v.nmemb++;
	}

	qsort(v.arr, v.nmemb, sizeof(*v.arr), cmp_double);


	size_t index = 0;

	// quantile
	for(size_t i = 0; i <= nquan; i++)
	{
		double perc = 1.0 * i / nquan;
		index = (size_t)(perc * v.nmemb);
		if(index >= v.nmemb) index = v.nmemb - 1;

		fprintf(stdout, "%g\t%g\n", perc, v.arr[index]);
	}

	// histogram
	double upper_bound = first_lower_bound + bin_width;
	size_t counter = 0;
	for(size_t i=0; i < v.nmemb; i++)
	{
		double value = v.arr[i];
		if(value<upper_bound) counter ++;
		else
		{
			fprintf(stderr, "%g\t%g\t%ld\n", upper_bound - bin_width, upper_bound, counter);
			counter = 0;
			upper_bound += bin_width;
			i--;
		}
	}
	free(v.arr);
	v.arr=NULL;
	free(linePtr);
	linePtr=NULL;

}
