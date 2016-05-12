VCFStatsAlive
=============

A utility that calculates statistics off a vcf stream, reporting at a given
interval in json

Note: Please use the stable branch, as the master branch contains untested changes.

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

If no vcf-file is specified, input is then read from stdin
```
