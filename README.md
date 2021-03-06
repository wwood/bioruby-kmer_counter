# bio-kmer_counter

[![Build Status](https://secure.travis-ci.org/wwood/bioruby-kmer_counter.png)](http://travis-ci.org/wwood/bioruby-kmer_counter)

bio-kmer_counter is a simple [biogem](http://biogem.info) for fingerprinting
nucleotide sequences by counting the occurences of particular kmers in the
sequence. The methodology is not new, for a reference see 
[Teeling et. al. 2004](http://www.biomedcentral.com/1471-2105/5/163). 
The default parameters are derived from the well explained methods section of
[Dick et. al. 2009](http://genomebiology.com/content/10/8/R85).

This methodology is quite different to that of other software that counts
kmer content with longer kmers, e.g. [khmer](https://github.com/ged-lab/khmer).
Here only small kmers are intended (e.g. 1-mer or 4-mer).

## Installation

After installing [Ruby](http://www.ruby-lang.org) itself, install the bio-kmer_counter rubygem:

```sh
gem install bio-kmer_counter
```

bio-kmer_counter is only tested on Linux, but probably works on OSX too. It might even work on Windows if
the progress bar is turned off. Maybe.

## Usage

The default parameters analyse a fasta file that contains one or more sequences in it for 4-mer (tetranucleotide)
content. By default, any sequence 
in the fasta file 2kb or longer is included at least once. Sequences are split up
into 5kb windows if they are that long, and each window is reported separately.
If the leftover bit at the end after any 5kb windows is 2kb or longer then this is also included.

By default, each 4 base window in the input sequence is included exactly once in the output file.
To account for the fact 
that the directions of sequences with respect to each other are presumed to be unknown (as is the
case for de-novo genome assembly), either the forward or reverse complement is included. Which one
(forward or reverse) depends on which one comes first alphabetically. So for instance if the window is ```CTTT```, then ```AAAG```
is used. Accounting for palindromic sequences like ```ATAT```, there are 136 of these lowest lexigraphical 4-mers.
So there are 136 columns in the output, plus one for the name of the window. Using only 1 is
actually slightly different than the method outlined in Dick et. al. 2009, but we
don't expect the results to differ.

Example usage, if you wish to fingerprint a fasta file ```my_nucleotide_sequences.fasta```:
```sh
kmer_counter.rb my_nucleotide_sequences.fasta >tetranucleotide_content.csv
```

The fingerprints are reported in percentages. Well, between 0 and 1, that is.
From there it is up to you how to use the fingerprints, sorry. For the full
gamut of options, see

```sh
kmer_counter.rb -h
```

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

