all: tests ibdtools

clean: 
	rm -rf ./gtest.dSYM ./.cache ./gtest  ./tmp*

tests: \
	src/gtest.cpp \
	include/chromosomes.hpp \
	include/common.hpp \
	include/genotypes.hpp \
	include/gmap.hpp \
	include/ibdcoverage.hpp \
	include/ibdfile.hpp \
	include/ibdmatrix.hpp \
	include/ibdmerger.hpp \
	include/ibdsorter.hpp \
	include/ibdspliter.hpp \
	include/ibdstat.hpp \
	include/metafile.hpp \
	include/positions.hpp \
	include/samples.hpp 

	${CXX} -std=c++17 -I src/cxxopts/include ${CXXFLAGS} -lz -lhts -lfmt -lgtest  -g  src/gtest.cpp  -o tests

ibdtools: \
	src/ibdtools.cpp \
	include/chromosomes.hpp \
	include/common.hpp \
	include/genotypes.hpp \
	include/gmap.hpp \
	include/ibdcoverage.hpp \
	include/ibdfile.hpp \
	include/ibdmatrix.hpp \
	include/ibdmerger.hpp \
	include/ibdsorter.hpp \
	include/ibdspliter.hpp \
	include/ibdstat.hpp \
	include/metafile.hpp \
	include/positions.hpp \
	include/samples.hpp 

	${CXX} -std=c++17 -I src/cxxopts/include ${CXXFLAGS} -lz -lhts -lfmt -O3 src/ibdtools.cpp -o ibdtools 
