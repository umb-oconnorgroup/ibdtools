FLAGS = -g -Wall -lm -lpthread

argtable3.o: argtable3.h argtable3.c
	gcc -c argtable3.c -g -Wall

vector.o: vector.h vector.c
	gcc -c vector.c -g -Wall

tpool.o: tpool.h tpool.c
	gcc -c tpool.c -g -Wall

ibdqc.o: ibdqc.c
	gcc -c ibdqc.c -g -Wall

ibdqc: ibdqc.o tpool.o vector.o argtable3.o
	gcc -g -o ibdqc ibdqc.o tpool.o vector.o argtable3.o $(FLAGS)

clean: 
	rm *.o ibdqc

