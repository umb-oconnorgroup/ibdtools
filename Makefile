INC_PATH = -I/usr/local/packages/htslib-1.11/include -I/usr/local/packages/tbb/include 
LIB_PATH = -L/usr/local/packages/htslib-1.11/lib  -L/usr/local/packages/tbb/lib/intel64/gcc4.4 
LIBS = -lz -lhts -ltbb
CC = g++
FLAGS = -std=c++17 -O3
FLAGS2 = -std=c++17 -g

build/ibdmerge : src/ibdmerge.cpp
	if [ ! -d build ]; then mkdir build; fi
	${CC} ${FLAGS} src/ibdmerge.cpp -o build/ibdmerge ${INC_PATH} ${LIB_PATH} ${LIBS}

cmp: build/ibdmerge
	if [ ! -d data ]; then mkdir data; fi
	if [ ! -f data/merge-ibd-segments.17Jan20.102.jar ]; then \
		wget https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar \
		-P data; \
	fi
	time build/ibdmerge test/example.ibd test/example.vcf.gz test/example.map -m 0.6 -d 1 \
		> test/merged.ibd
	./test/cmp2browning.sh

