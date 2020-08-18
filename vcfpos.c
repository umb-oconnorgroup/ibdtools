#include<stdio.h>
#include<stdlib.h>

int main()
{
	char * buff = calloc(1000, 1);;
	size_t buff_size = 999;
	char *p=NULL, *p2=NULL;;

	while((getline(&buff, &buff_size, stdin) > 0))
	{
		// printf("%ld\n", buff_size);
		if(*buff == '#') continue;
		else{
			p = buff;
			while(*p != '\t')p++;
			p++;

			p2 = p;
			while(*p2 != '\t')p2++;
			p2++;
			*p2 = '\0';

			fprintf(stdout, "%s\n", p);
		}
	}

	if(buff != NULL){
		free(buff);
		buff=NULL;
	}
	p = p2 = NULL;
	return 0;
}
