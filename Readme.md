# `ibdmerge`

`ibdmerge` merges IBD segments following Browning's tool `merge-ibd-segments.jar`.
This is written to handle large datasets.

## compilation

1. Adjust the paths for libraries: `hts`, `zlib` and `tbb` in the `Makefile`; 
2. Add these lib paths to the `LD_LIBRARY_PATH` environment variable if not found in
   default paths.
3. Run `make` at the root directory to compile. The binaries will be under `build/`.

## Usage:	

```
Usage: ibdmerge [OPTION...] IBD_FILE VCF_FILE MAP_FILE

 Required arguments:
  IBD_FILE                   IBD files with format `sn1:sn2     start   end...`,
                             sorted by sort -k1,1 -k2,2n -k3,3n
  MAP_FILE                   MAP file used to call IBD
  VCF_FILE                   VCF file used to call IBD. The VCF file is
                             expected to only have biallelic sites

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

