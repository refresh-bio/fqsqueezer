# FQSqueezer - FASTQ compressor based on k-mer statistics

FQSqueezer is an experimental high-end compression of short-read FASTQ files. 
The main goal of the tool is to offer the best possible compression ratio with running times allowing to run it even for WGS human datasets.

FQSqueezer usually offers compression ratios tens of percent better than given by the state-of-the-art tools, like FaStore, Minicom, Spring. The running times are, however, significantly longer.

## Quick start
```
# Clone and make
git clone https://github.com/refresh-bio/fqsqueezer.git
cd fqsqueezer && make -j

# Compress a single FASTQ file using 16 threads with read reorganization
./fqs-1.1 e -s -t 16 -out SRR105788_1.fqs SRR105788_1.fastq

# Compress a single FASTQ file using 16 threads without read reorganization
./fqs-1.1 e -s -om o -t 16 -out SRR105788_1.fqs SRR105788_1.fastq

# Decompress a single FASTQ file
./fqs-1.1 d -out SRR105788_1.fastq SRR105788_1.fqs

# Compress paired-end FASTQ files using 16 threads with read reorganization
./fqs-1.1 e -p -t 16 -out SRR105788.fqs SRR105788_1.fastq SRR105788_2.fastq

# Compress paired-end FASTQ files using 16 threads with read reorganization
./fqs-1.1 e -p -t 16 -om o -out SRR105788.fqs SRR105788_1.fastq SRR105788_2.fastq

# Decompress paired-end FASTQ files
./fqs-1.1 d -out SRR105788_1.fastq -out2 SRR105788_2.fastq SRR105788.fqs
```

## Installation
FQSqueezer can be downloaded from this repo and compiled. 
The supported PS are:
* Windows: Visual Studio solution provided,
* Linux: make project (G++ 9.0 or newer required).

## Version history 
* 1.1 (16 June 2020)
  * bugfix release
* 1.0 (24 February 2019)
  * first public release
 
## Usage
### Compression
```
fqs-1.1 e [compression-options] <input.fastq>
fqs-1.1 e [compression-options] <input1.fastq> <input2.fastq>
fqs-1.1 e [compression-options] @<file_list>
```

Compression options:
* `-s` &ndash; single-end data
* `-p` &ndash; paired-end data
* `-t <num>` &ndash; number of threads (default: `1`)
* `-gs <num>` &ndash; approx. length of genome in Mbp (default: `3100`)
* `-tmp <path>` &ndash; path to temporary files (default: `./fqs_tmp_`)
* `-out <path>` &ndash; path to output file (default: `output.fqs`)
* `-om <s|o>` &ndash; order of reads:
  * `s` &ndash; sorted (default) 
  * `o` &ndash; original
* `-qm <o|8|4|2|n>` &ndash; quality mode:
  * `o` &ndash; original
  * `8` &ndash; Illumina 8-lev. (default)
  * `4` &ndash; Illumina 4-lev.
  * `2` &ndash; binary thr
  * `n` &ndash; none
* `-qt <num>` &ndash; threshold of quality for "trusted" base (default: `20`)
* `-im <o|i|n>` &ndash; id mode:
  * `o` &ndash; original
  * `i` &ndash; instrument only (default)
  * `n` &ndash; none
* `-v <num>` &ndash; verbosity from 0 to 2 (default: `1`)

### Decompression:
```
fqs-1.1 d [decompression-options] <input.fastq>
```

Decompression options:
* `-out <path>` &ndash; path to 1st (or only in SE mode) output file (default: `output.fqs`)
* `-out2 <path>` &ndash; path to 2nd output file in PE mode (default: `output2.fqs`)

## Citing
[Deorowicz, S. (2020) FQSqueezer: k-mer-based compression of sequencing data, Scientific Reports, 10:578.](https://www.nature.com/articles/s41598-020-57452-6)
