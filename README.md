[TOC]

LocARNA: Alignment of RNAs
==========================

LocARNA is a collection of alignment tools for the structural analysis
of RNA. Given a set of RNA sequences, LocARNA simultaneously aligns
and predicts common structures for your RNAs. In this way, LocARNA
performs Sankoff-like alignment and is in particular suited for
analyzing sets of related RNAs without known common structure.

LocARNA distinguishes itself from many other Sankoff-style multiple
alignment programs by its performance and low memory complexity, high
accuracy, and richness of features. As unique features, it offers
structure-local alignment, flexible structure and anchor constraints,
and provides efficient computation of reliabilities in
sequence-structure alignment. The package offers a robust core of
features and is used as experimental platform for the incorporation of
new features in RNA sequence-structure alignment.


At its core, the package offers global and local multiple alignment of
RNAs.

Multiple alignment can be performed in one of several different ways:

* progressive alignment using sequence-structure alignment of profiles

* progressive alignment after consistency transformation using
T-Coffee

* progressive alignment using probabilistic consistency transformation
  and sequence-structure profile alignments, optionally followed by
  iterative refinement.


Besides of global alignment, LocARNA supports two kinds of
locality. Local alignment as it is known from sequence alignment,
identifies and aligns the best matching subsequences. This form of
locality is called sequence local to distinguish it from structural
locality. When performing structure local alignment, LocARNA
identifies and aligns the best matching substructures in the RNAs. The
sequences of those substructures can be discontinuous on the sequence
level, but remain connected via structural bonds.

Alignment Reliabilities (LocARNA-P). In this special, probabilistic
mode of operation LocARNA supports the efficient computation of match
probabilities, probabilistic consistency transformation for more
accurate multiple alignment, and generates reliability profiles of
multiple alignments.


------------
Installation
------------

The software can be installed on recent Linux or MacOS X systems. Windows
is untested but should be supported via Linux extensions.

### Installation from Conda package (recommended)

On Mac/Linux, LocARNA is installed mosteasily via Conda from a pre-compiled
package. For this purpose, install Conda and run

```
conda install -c bioconda locarna
```

from the command line.


### Alternative installation from source

Installing from source requires a C++ build system: typically, g++ or
clang and autotools. Moreover, it depends on the Vienna RNA package.

#### Installation from source distribution

Obtain the tar.gz source distribution, e.g. from Github

[https://github.com/s-will/LocARNA/releases](https://github.com/s-will/LocARNA/releases)

Then, build and install like

```
tar xzf locarna-xxx.tar.gz
cd locarna-xxx
./configure --prefix=/usr/local
make
make install
```

Is Vienna RNA installed in a non-standard location, this has to be
specified by configure option ```--with-vrna=path-to-vrna``.

Installing from source furthermore allows testing via

```
make test
```

and building documentation locally by

```
make doxygen-doc
```

Building documentation requires additional tools: doxygen, pod2markdown and pandoc.



#### Installation from the git repository

Installing from repository is possible after cloning and setting up the
autotools suite. This is most easily achieved by running

```
autoreconf -i
```

in the cloned repository. Then, the installation essentially works like
installing from source distribution. Note that, we will however require
additional tools to build the documentation: help2man, pod2man.



-----
Usage
-----

For instructions on the use of the tools, please see the documentation / man pages of
the single tools

* [mlocarna](md_src_Utils_mlocarna.html) --- for multiple alignment of
  RNAs.  This program supports most of the functionality of the package via
  a high level interface.

* [locarna](md__build_Doc_man_locarna.html) --- for pairwise alignment

* [locarna_p](md__build_Doc_man_locarna_p.html) --- for pairwise computation of alignment partition function
  and (sequence and structure) match probabilities

* [sparse](md__build_Doc_man_sparse.html) --- for structurally stronger sparsified pairwise alignment


For additional functionality and special purposes, see

* [exparna_p](md__build_Doc_man_exparna_p.html) --- for generating exact matches from the ensembles of two RNAs

* [locarnate](md_src_Utils_locarnate.html) --- for multiple alignment of
  RNAs via T-Coffee. This script offers multiple alignment of RNAs that is
  performed by sequence-structurally aligning all pairs of RNAs and then
  using T-Coffee to construct a common multiple alignment out of all
  pairwise ones.


-------------------------------------------
Online information, download and web server
-------------------------------------------

Further information is provided online at

   http://www.bioinf.uni-freiburg.de/Software/LocARNA/

The most recent versions of the package will be made available through
this page.


The core functionality of the package is accessible through a web
interface at

   http://rna.informatik.uni-freiburg.de


-------
Contact
-------

Main author and contact: Sebastian Will sebastian.will (at) polytechnique.edu
