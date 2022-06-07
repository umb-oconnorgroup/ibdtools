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
most efficient, but has been used in real research projects. As projects progress,
we will add more functions into this program.

- `ibdtools encode`: encode the IBD file, VCF file and plink map file into
  binary format for better/quicker IO
- `ibdtools split`: remove IBD segments overlapping with given regions or
  calculated regions with low SNP density. This can be useful for removing false
  positive IBD segments in regions of low SNP density, or correcting
  selection-induced bias.
- `ibdtools sort`: sort large IBD files by implementing an external sorting
  algorithm. It must be run before calling `ibdtools merge`.
- `ibdtools merge`: flatten haplotype-pair IBD segments into individual-pair
  IBD segments by merging all IBD segments shared by a pair of individuals when
  they overlap or are separated by a short gap with only a few discordant
  sites. This sub-command reimplements the algorithm in Dr. Browning's
  `merge-ibd-segments` tool for stability and consistency purposes.
- `ibdtools matrix`: allow aggregating IBD segments into chromosome-wide and
  genome-wide total IBD and output total IBD matrix for downstream analysis.
  During the aggregation process, IBD can be filtered at different levels
  including the subpopulation, IBD segment length, and total IBD length.
- `ibdtools snpdens`: calculate SNP density across chromosome.
- `ibdtools coverage`: calculate IBD coverage across chromosome.
- `ibdtools view`: search and print IBD segments belonging to a specified pair of
  individuals.
- `ibdtools decode`: convert the processed, encoded IBD back to text format for
  readability and downstream analysis.
- `ibdtools stat`: currently, only calculate the IBD length distribution.

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
