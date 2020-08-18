#include<stdio.h>

int main()
{
	char buff[1000];
	size_t buff_size = 1000;
	char *linePtr = buff, *p=NULL, *p2=NULL;;

	while((getline(&linePtr, &buff_size, stdin) > 0))
	{
		if(*linePtr == '#') continue;
		else{
			p = linePtr;
			while(*p != '\t')p++;
			p++;

			p2 = p;
			while(*p2 != '\t')p2++;
			p2++;
			*p2 = '\0';

			fprintf(stdout, "%s\n", p);
		}
	}
	return 0;
}
