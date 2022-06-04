all: tests ibdtools

clean: 
	rm -rf ./gtest.dSYM ./.cache ./gtest  ./tmp*

tests: \
	gtest.cpp \
	../include/chromosomes.hpp \
	../include/common.hpp \
	../include/genotypes.hpp \
	../include/gmap.hpp \
	../include/ibdcoverage.hpp \
	../include/ibdfile.hpp \
	../include/ibdmatrix.hpp \
	../include/ibdmerger.hpp \
	../include/ibdsorter.hpp \
	../include/ibdspliter.hpp \
	../include/ibdstat.hpp \
	../include/metafile.hpp \
	../include/positions.hpp \
	../include/samples.hpp 

	${CXX} -std=c++17 -I ../include -I cxxopts/include ${CXXFLAGS} -lz -lhts -lfmt -lgtest  -g  gtest.cpp  -o tests

ibdtools: \
	ibdtools.cpp \
	../include/chromosomes.hpp \
	../include/common.hpp \
	../include/genotypes.hpp \
	../include/gmap.hpp \
	../include/ibdcoverage.hpp \
	../include/ibdfile.hpp \
	../include/ibdmatrix.hpp \
	../include/ibdmerger.hpp \
	../include/ibdsorter.hpp \
	../include/ibdspliter.hpp \
	../include/ibdstat.hpp \
	../include/metafile.hpp \
	../include/positions.hpp \
	../include/samples.hpp 

	${CXX} -std=c++17 -I ../include -I cxxopts/include ${CXXFLAGS} -lz -lhts -lfmt -O3 ibdtools.cpp -o ibdtools 
