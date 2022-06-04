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

	${CXX} -std=c++17 \
		-I src/cxxopts/include ${CXXFLAGS} \
		-lz -lhts -lfmt -lgtest  -g  \
		src/gtest.cpp  -o tests

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

	${CXX} -std=c++17 \
		-I src/cxxopts/include ${CXXFLAGS} \
		-lz -lhts -lfmt -O3 \
		src/ibdtools.cpp -o ibdtools 

ibdtools.AppImage: ibdtools

	if [ ! -e ./linuxdeploy-x86_64.AppImage ]; then \
		wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage; \
		chmod +x ./linuxdeploy-x86_64.AppImage; \
	fi
	./linuxdeploy-x86_64.AppImage \
		-e ibdtools \
		-l $$(ldd ibdtools | grep libstdc++ | sed -e 's/^\s*//' | awk '{print $$3}') \
		-l $$(ldd ibdtools | grep libhts    | sed -e 's/^\s*//' | awk '{print $$3}') \
		-l $$(ldd ibdtools | grep libfmt    | sed -e 's/^\s*//' | awk '{print $$3}') \
		--appdir ibdtools.AppDir \
		-d appimage/ibdtools.desktop \
		-i appimage/icon_256x256.png \
		--output appimage
	rm -rf ibdtools.AppDir
	mv ibdtools*AppImage ibdtools.AppImage

