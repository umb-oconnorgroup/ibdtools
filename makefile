FLAGS = -g -Wall -lm -lpthread

vector.o: vector.h vector.c
	gcc -c vector.c -g -Wall

tpool.o: tpool.h tpool.c
	gcc -c tpool.c -g -Wall

ibdqc.o: ibdqc.c
	gcc -c ibdqc.c -g -Wall

ibdqc: ibdqc.o tpool.o vector.o
	gcc -g -o ibdqc ibdqc.o tpool.o vector.o $(FLAGS)

clean: 
	rm *.o ibdqc

