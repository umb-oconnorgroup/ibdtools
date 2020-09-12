#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>

int main(int argc, char * argv[])
{
	char * lineBuff =NULL, *tok=NULL;
	size_t lineBuffSize=0;

	if(argc != 2)
	{
		fprintf(stderr, "\n\n Usage: cat xxx.ibd | ./ibdhap2dip <nsam> > xxx_dip.ibd\n\n"
				"\n input ibd should have at least 9 columns, 8th col=score, 9th col= ibd length\n");
		exit(0);
	}
	int nsam = atoi(argv[1]);

	while(getline(&lineBuff, &lineBuffSize, stdin)>0 && *lineBuff != '\n')
	{
		tok = strtok(lineBuff, " \t\n");
		int id1 = atoi(tok);
		tok = strtok(NULL, " \t\n");
		int hap1 = atoi(tok);
		tok = strtok(NULL, " \t\n");
		int id2 = atoi(tok);
		tok = strtok(NULL, " \t\n");
		int hap2 = atoi(tok);
		tok = strtok(NULL, " \t\n");
		int chr = atoi(tok);
		tok = strtok(NULL, " \t\n");
		size_t start = strtol(tok, NULL, 10);
		tok = strtok(NULL, " \t\n");
		size_t end = strtol(tok, NULL, 10);
		tok = strtok(NULL, " \t\n");
		assert(tok != NULL);
		double score = strtod(tok, NULL);
		tok = strtok(NULL, " \t\n");
		assert(tok != NULL);
		double cM = strtod(tok, NULL);

		if(id1 >= nsam) {	
			id1 -= nsam;
			hap1 = 2;
		}
		if(id2 >= nsam) {
			id2 -= nsam;
			hap2 = 2;
		}

		if(id1>id2)
		{
			int temp = id1;
			id1 = id2;
			id2 = temp;
		}

		fprintf(stdout, "\%d\t%d\t%d\t%d\t%d\t%ld\t%ld\t%g\t%g\t\n",
				id1, hap1, id2, hap2, chr, start, end, score, cM);
	}

	return 0;
}
