#! /bin/bash

set -e

# use msprime to simulate vcf and call ibd using hapibd. The output is used as input for ibdtools
python3 ./simulate_data.py

# use ibdtools to encode, split, deocde, sort, merge ibd and calculate chromosome-level total ibd matrix
for i in {1..5}; do 
	../build/ibdtools encode -i $i.ibd.gz -v $i.vcf.gz -g $i.map -c $i -o $i.eibd -m $i.meta -M 1
	../build/ibdtools split -i $i.eibd -m $i.meta -W 2 -S 10 -o $i.split -M 1
	../build/ibdtools decode -i $i.split1 -m $i.meta -o $i.ibd.gz -M 1
	../build/ibdtools sort -i $i.split1 -o $i.sibd -M 1
	../build/ibdtools merge -i $i.sibd -m $i.meta -o $i.mibd -M 1
	../build/ibdtools matrix -i $i.mibd -m $i.meta -M 1 -T 2.5 -L 5 -o $i.mat
done

# collect all chromosome-wide total matrices and join them with comma
MATRICES=`\ls {1..5}.mat*.mat | tr '\n' ','| sed s/,$//` 
../build/ibdtools matrix -x $MATRICES -m 1.meta -T 2.5 -L 5 -o gw.mat -M 3
