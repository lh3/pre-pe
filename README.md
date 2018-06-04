Pre-pe is a set of tools to preprocess paired-end reads. They trim generic
sequencing adapters, clip experiement specific adapters, identify inline
barcodes and merge overlapping ends. Each tool in pre-pe has very similar
functionality but for different data types (see below). The general command
line to run pre-pe looks like:

```sh
seqtk mergepe read1.fq.gz read2.fq.gz | pre-adna - | gzip > pe-se-mixed.fq.gz
```
Here, the [seqtk][seqtk] command line generates an interleaved FASTQ stream.
You can skip this step if your FASTQs are already interleaved. By default,
paired-end and merged single-end reads are merged into a single stream. You can
use [bwa-mem][bwa] to directly map the output with

```sh
pre-meta interleaved.fq.gz | bwa mem -p ref.fa - | gzip > output.sam.gz
```
Most other short-read mappers don't have this functionality.

Pre-pe consists of the following tools:

* pre-lianti for single-cell whole-genome sequencing data produced with the
  [LIANTI protocol][lianti]. It is the first tool in this series and also
  available from my [lianti repo][lianti-repo].

* pre-adna for ancient data produced with the [Reich Lab protocol][udg]. This
  program is also available from the [adna repo][adna]. It has been used to
  process hundreds of ancient full genomes from the Reich lab.

* pre-dip-c for single-cell Hi-C data produced with the Dip-C protocol. It was
  modified from pre-lianti by [Longzhi Tan][longzhi].

* pre-meta for single-cell genomic data produced with the META protocol. It is
  similar to pre-dip-c except that it additionally checks ambiguous end merges
  to greatly reduce artifactual deletions.

Of these tools pre-adna may be of general interest.

[seqtk]: https://github.com/lh3/seqtk
[bwa]: https://github.com/lh3/bwa
[udg]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4275898/
[adna]: https://github.com/DReichLab/adna
[lianti]: https://www.ncbi.nlm.nih.gov/pubmed/28408603
[lianti-repo]: https://github.com/lh3/lianti
[longzhi]: https://github.com/tanlongzhi
