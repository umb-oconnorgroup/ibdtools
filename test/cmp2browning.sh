#! /bin/bash

CHR=10
# 1 reformat for input and run hap-ibd
time cat test/example.ibd | tr ':' '\t' \
	| awk -v OFS='\t' -v chr=$CHR '{print $1, 1, $2, 2, chr, $3, $4, $5}' \
	\
	\
	| java -jar data/merge-ibd-segments.17Jan20.102.jar \
	test/example.vcf.gz \
	test/example.map \
	0.6 1 \
	> test/merged_browning.ibd

cut -f1-4  test/merged.ibd | sort -k1,1 -k2,2 -k3,3n -k4,4n > test/cmp_1.txt
cut -f 1,3,6,7 test/merged_browning.ibd | sort -k1,1 -k2,2 -k3,3n -k4,4n > test/cmp_2.txt

diff test/cmp_1.txt test/cmp_2.txt

if [ $? == 0 ]; then
	echo 'ibdmerge and merge-ibd-segments.jar generated the SAME result!'
fi
rm test/cmp_1.txt test/cmp_2.txt

