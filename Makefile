CC = ${HOME}/miniconda3/bin/x86_64-conda_cos6-linux-gnu-g++
LIBS = -lz -lhts -ltbb
FLAGS = -std=c++17 -O3
FLAGS2 = -std=c++17 -g

# check default include path by run 
# xxx_gcc -x c++ -v -E /dev/null
INC_PATH = -I${HOME}/miniconda3/include/ \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/../../../../x86_64-conda-linux-gnu/include/c++/9.3.0 \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/../../../../x86_64-conda-linux-gnu/include/c++/9.3.0/x86_64-conda-linux-gnu \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/../../../../x86_64-conda-linux-gnu/include/c++/9.3.0/backward \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/include \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/include-fixed \
 -I${HOME}/miniconda3/bin/../lib/gcc/x86_64-conda-linux-gnu/9.3.0/../../../../x86_64-conda-linux-gnu/include \
 -I${HOME}/miniconda3/bin/../x86_64-conda-linux-gnu/sysroot/usr/include

LIB_PATH = -L${HOME}/miniconda3/lib



all: build/ibdsort build/ibdmerge

build/ibdsort : src/ibdsort.cpp
	if [ ! -d build ]; then mkdir build; fi
	${CC} ${FLAGS} src/ibdsort.cpp -o build/ibdsort ${INC_PATH} ${LIB_PATH} ${LIBS}

build/ibdmerge : src/ibdmerge.cpp
	if [ ! -d build ]; then mkdir build; fi
	${CC} ${FLAGS} src/ibdmerge.cpp -o build/ibdmerge ${INC_PATH} ${LIB_PATH} ${LIBS}

build/ibdtotal : src/ibdtotal.cpp
	if [ ! -d build ]; then mkdir build; fi
	${CC} ${FLAGS} src/ibdtotal.cpp -o build/ibdtotal ${INC_PATH} ${LIB_PATH} ${LIBS}

test_ibdsort: build/ibdsort
	build/ibdsort test/example_shuffled.ibd -s test/example_sample_name.txt -m 0.000016000 -k 3  > ./test/example_sort.ibd
	awk -v OFS='\t' '{if($$1 > $$3){tmp=$$3; $$3=$$1; $$1=tmp;} print $$1, $$3, $$6, $$7}' ./test/example_shuffled.ibd | sort -k1,1 -k2,2 -k3,3n -k4,4n > test/example_sort_cmd.ibd
	diff test/example_sort.ibd test/example_sort_cmd.ibd
	@if [ $$? == 0 ]; then echo; echo; echo ibdsort result is the SAME as sort command; echo; echo; fi


test_ibdmerge: build/ibdmerge
	if [ ! -d data ]; then mkdir data; fi
	if [ ! -f data/merge-ibd-segments.17Jan20.102.jar ]; then \
		wget https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar \
		-P data; \
	fi
	time build/ibdmerge test/example.ibd test/example.vcf.gz test/example.map -m 0.6 -d 1 \
		> test/merged.ibd
	./test/cmp2browning.sh

test_ibdtotal: build/ibdtotal
	build/ibdtotal test/example_sample_name.txt test/total.out -i test/example.ibd 
	build/ibdtotal test/example_sample_name.txt test/total.out2 -i test/example.ibd \
		-m test/total.out
	test/test_ibdtotal.py
