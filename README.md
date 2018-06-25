VCFStatsAlive
=============

A utility that calculates statistics off a vcf stream, reporting at a given
interval in json

Note: Please use the stable branch, as the master branch contains untested changes.
Note: autoconf/autoheader must be installed in order to build htslib dependency
Note: On linux systems, use 'sudo make' rather than 'make'. Also export LD_LIBRARY_PATH="/usr/local/lib" before
running so that htslib can be found for dynamic loading.

Usage
=====

```
vcfstatsalive [options] [vcf-file]

Options:
  -u	updateRate [default=1000]		The number of reads vcfstatsalive needs to process before producing another statistics update
  -f	firstUpdateRate [default=0]		The number of reads vcfstatsalive needs to process before producing the first statistics update. Useful to increase app responsiveness
  -q	qualHistLowerVal [default=1]	The lower value of invalid QUAL value. Any QUAL value less than this will not be counted towards quality histogram
  -Q	qualHistUpperVal [default=200]	The upper value of invalid QUAL value. Any QUAL value greater than this will not be counted towards quality histogram
  -l	logScaleAF [default=false]	    When specified, allele frequency histogram will be in log scale
  -b	batch [default=false]	    When specified, the statistics will only be outputed a single time at the end of the analysis.
<<<<<<< HEAD
  -S	subsetsamples [default=false]	    When specified for a multisample vcf the stats will be calculated for each sample independently and the output will be in the form `sampleName.json` with 1 file for each sample.
=======
>>>>>>> master

If no vcf-file is specified, input is then read from stdin
```
