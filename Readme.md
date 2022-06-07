# Introduction

An identical-by-descent (IBD) segment is an identical-by-state (IBS) segment
shared by a pair of haplotypes/individuals that is inherited from a common ancestor
without recombination. The concept is important in understanding the genetic
relatedness, population structure and effective population size in recent
generations. As different individual pairs may share IBD segments of different
lengths over different genomic regions, the number of IBD segments increase
quadraticlly with increasing sample size, and proportionally with genome size. 
High recombination rate as well as genotyping/phasing errors can break down long
IBD segments into multiple shorter segments, making the scale even larger.

Recently, several IBD calling tools, such as hap-IBD, iLASH, and TPBWT, have
been designed to efficiently detect IBD segments using phased genotype data.
However, fewer tools were designed to handle this large-scale IBD data. The
convenience of R and python in handling the IBD result (table/data frame) comes
with a cost of memory and a possibly slower computation speed. Depending
on the computing environment, high requirement of memory and dynamically
reallocating memory might lead to instability of these programs.

`ibdtools` is an attempt to improve the efficiency of handling large scale IBD
segments information by:
- better encoding the IBD segments
- stabler and tighter memory management 
- implementing some commonly-used algorithms in c/c++

Currently, the functionalities and algorithms are very basic and probably not the
most efficient, but has been to real research projects. As projects progress,
we will add more functions into this program.

# Clone the repo and the submodules
```sh
git clone git@github.com:gbinux/ibdtools.git --recurse-submodules
cd ibdtools
```
# Dependencies
1. `htslib`
2. `fmt` 
3. `gtest`

They can be install using conda
```sh
conda env create -f env.yml
```
# Compile ibdtools

```sh
conda activate ibdtools
cd src
make ibdtools
```

# Simulated data and example

```sh
cd example

# simulating data
./simulate_data.py

# running ibdtools on the simulated input data
example/example.sh
```

# Documentation

1. List all available subcommand:
```sh
./ibdtools 
```

2. Show help message/documentation for a subcommand
```sh
./ibdtools [subcommand] -h
```
