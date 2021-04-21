# ibdtools

## Compilation

1. Dependencies: The code relies on zlib, htslib and tbb libraries and  a version gcc
   compilers that supports C++ 17 standard. They can be installed by conda.

```sh
conda install -c conda-forge gxx_linux-64 htslib tbb-devel 
```

2. Adjust the paths for libraries: `hts`, `zlib` and `tbb` in the `Makefile`; 

3. Add these lib paths to the `LD_LIBRARY_PATH` environment variable if not found in
   default paths.

4. Run `make` at the project root directory to compile. The binaries will be under
   `build/`.


## `ibdsort`

```
Usage: ibdsort [OPTION...] IBD_FILE
ibdsort prepares formatted and sorted IBD data for ibdmerge. It uses external merge sort
algorithm and is written to handle large datasets. 

 Required arguments:
  IBD_FILE                   Input IBD file from hap_ibd for a single
                             chromosome

 Optional arguments:
  -k, --k_ways=integer       Number of ways for merging. Default 8
  -K, --keep_tmp_file        Keep all tmp file for debuging
  -m, --memory=float         Estimated maximum memory in Gb to use. Default 1
  -o, --out=filename         IBD output file after formatting and sorting,
                             default stdout
  -p, --tmp_file_prefix=string   The tmp file prefix. The file name will be
                             {tmp_file_prefix}_{xx}. Default is using the
                             IBD_FILE.
  -s, --sample=filename      Sample file. If provided, IBD records will be
                             sorted according to the alphabetical order of
                             sample name; if not provodied, IBD samples will be
                             sorted according to the order the program first
                             encounters the sample name in the IBD records
                             Default: NULL
  -v, --verbose              Show processing updates, default false

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

```

## `ibdmerge`

```
Usage: ibdmerge [OPTION...] IBD_FILE VCF_FILE MAP_FILE
ibdmerge merges IBD segments following Browning's tool merge-ibd-segments.jar.
This is written to handle large datasets.

 Required arguments:
  IBD_FILE                   IBD files with format `sn1:sn2	start	end...`,
                             sorted by sort -k1,1 -k2,2n -k3,3n
  MAP_FILE                   MAP file (space delimited) used to call IBD
  VCF_FILE                   VCF file (tab delimited) used to call IBD. The VCF
                             file is expected to only have biallelic sites

 Optional arguments:
  -d, --discord=int          max number of genotypes in IBD gap that are
                             inconsistent with IBD, default 1
  -m, --max_cM=float         max length of gap between IBD segments (cM),
                             default 0.6
  -o, --out=filename         IBD output file after merging, default stdout
  -v, --verbose              Show processing updates, default false

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```
