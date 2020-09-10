#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main()
{
    size_t size = 10000 * 10000;
    char *buff = calloc(size, sizeof(char));
    char *buff2 = calloc(size, sizeof(char));
    char *p = buff;
    char *ptemp;
    long col = 0, count;
    while (getline(&buff, &size, stdin) >= 0) {
        // header ##
        if (strncmp(buff, "##", 2) == 0)
            printf("%s", buff);
        // header #CHROM
        else if (strncmp(buff, "#C", 2) == 0) {
            for (char *p1 = buff; p1 < buff + strlen(buff); ++p1)
                if (*p1 == '\t')
                    ++col;
            ++col;
            printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for (int i = 0; i < col - 9; i++)
                printf("\t%d\t%d", 2 * i, 2 * i + 1);
            printf("\n");
        } else {
            for (char *p1 = buff, *p2 = buff2; p1 < buff + strlen(buff);) {
                if ((*p1 == '0' || *p1 == '1') && *(p1 + 1) == '|') {
                    p2[0] = p1[0];
                    p2[1] = '|';
                    p2[2] = p1[0];
                    p2[3] = '\t';
                    p2[4] = p1[2];
                    p2[5] = '|';
                    p2[6] = p1[2];
                    p2 += 7;
                    p1 += 3;
                } else {
                    *p2 = *p1;
                    ++p2;
                    ++p1;
                }
                *p2 = '\0';
            }
            printf("%s", buff2);
        }
    }
    return 0;
}
