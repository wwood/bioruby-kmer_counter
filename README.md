# bio-kmer_counter

[![Build Status](https://secure.travis-ci.org/wwood/bioruby-kmer_counter.png)](http://travis-ci.org/wwood/bioruby-kmer_counter)

bio-kmer_counter is a simple [biogem](http://biogem.info) for fingerprinting
nucleotide sequences by counting the occurences of particular kmers in the
sequence. The methodology is not new, for references see [Teeling et. al. 2004](http://www.biomedcentral.com/1471-2105/5/163). The default parameters are derived from the methods section of [Dick et. al. 2009](http://genomebiology.com/content/10/8/R85).

This methodology is quite different to that of other software that counts
kmer content with longer kmers, e.g. [khmer](https://github.com/ged-lab/khmer).
Here only small kmers are intended (e.g. 1mer or 4mer).

Note: this software is under active development!

## Installation

```sh
    gem install bio-kmer_counter
```

## Usage

To analyse a fasta file (that contains one or more sequences in it) for 4-mer (tetranucleotide)
content, reporting the fingerprint of 5kb windows in each sequence separately,
plus the leftover part if it is longer than 2kb:

```sh
    kmer_counter.rb <fasta_file> >tetranucleotide_content.csv
```

The fingerprints are reported in percentages. Well, between 0 and 1, that is.
From there it is up to you how to use the fingerprints, sorry.

## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/wwood/bioruby-kmer_counter

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

This software is currently unpublished, so please just cite the homepage (thanks!).

Please also cite the tools upon which it is based, one of:
  
* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

## Copyright

Copyright (c) 2012 Ben J Woodcroft. See LICENSE.txt for further details.

