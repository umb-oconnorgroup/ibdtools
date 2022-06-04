all: compile
	./gtest

debug: compile_debug
	lldb ./gtest

clean: 
	rm -rf ./gtest.dSYM ./.cache ./gtest  ./tmp*



compile_debug: 
	@g++ -std=c++20 \
		-isystem  /usr/local/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
		-isystem  /Library/Developer/CommandLineTools/usr/lib/clang/12.0.5/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include \
		-isystem  /Library/Developer/CommandLineTools/usr/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks \
		-L /usr/local/lib \
		-lz -lhts -lfmt -ltbb -lgtest \
		-g \
		gtest.cpp \
		-o gtest 
	

compile: 
	@g++ -std=c++20 \
		-isystem  /usr/local/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
		-isystem  /Library/Developer/CommandLineTools/usr/lib/clang/12.0.5/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include \
		-isystem  /Library/Developer/CommandLineTools/usr/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks \
		-L /usr/local/lib \
		-lz -lhts -lfmt -ltbb -lgtest \
		-L /usr/local/lib \
		-lz -lhts -lfmt -ltbb  \
		-O3 \
		 gtest.cpp \
		-o gtest

ibdtools: ibdtools.cpp
	@g++ -std=c++20 \
		-isystem  /usr/local/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
		-isystem  /Library/Developer/CommandLineTools/usr/lib/clang/12.0.5/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include \
		-isystem  /Library/Developer/CommandLineTools/usr/include \
		-isystem  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks \
		-L /usr/local/lib \
		-lz -lhts -lfmt -ltbb -lgtest \
		-L /usr/local/lib \
		-lz -lhts -lfmt -ltbb -lboost_program_options  \
		-O3 \
		 ibdtools.cpp \
		-o ibdtools 

test_ibdtools:
	./ibdtools encode \
		-c 10 \
		-g ../data/example.map \
		-v ../data/example.bcf.gz \
		-i ../data/example.ibd.gz \
		-o tmp.10.ibd \
		-m tmp.10.meta

