# ChangeLog

## 2.0.1
Maintenance updates

* Update the Catch testing framework to version 2.13.10
* Replace deprecated use of bind2nd to support recent C++ compilers


## 2.0.0

Version 2 of LocARNA introduces several improvements that make LocARNA
faster and easier to use. Most prominently, it features a novel sequence
envelope heuristic in all major alignment modes, which significantly
improves over less sophisticated banding heuristics (in terms of speed,
robustness, /and/ flexibility). Most directly notable for users, it changes
the default parameters to more appropriate settings. In addition to
various detailed improvements of the LocARNA tools, there are significant
changes under the hood: the C++-code of the LocARNA library was simplified
and modernized making use of C++-14 features. These changes make the
library easier to use, maintain and extend.

The development of version 2 extended over numerous 'release candidates'
that not only fixed bugs and stabilized features, but - despite their RC
designation - also introduced some new functionality e.g. (for
clustering and benchmarking). Comprehensive change
log is therefore provided for the sinlge release candidate releases.

---

## 2.0.0RC10

### Documentation and examples

* Improve readme
* Update doxygen documentation
* Fix: add conda-forge channel in conda install
* Add CLUSTAL-header in realign example
* Use full names in tRNA_5 example
* use doxygen awesome style
* use REAMDE.md as main page in docu
* Reformat changelog and include in docu
* include man pages from perl scripts and C++ programs
* Reformat help output of C++ programs
* Activate doxygen pdf docu
* Reformatting in mlocarna pod
* Update email address
* Remove build status

### Improvements of functionality / usability

* Make mlocarna less chatty
* mlocarna: Catch empty names in fasta input + modify output in case of duplicate names
* mlocarna: Catch missing clustal header in realignment input
* turn off warning about unsupported nucleotides
* Add benchmark plot script to compare multiple experiments
* Add alignments target in benchmark Makefile

### Fixes

* Keep input file in realignment mode as input.aln
* Fix pp2dot
* Fix reliability output (aligned offsets/name widths)
* Fix mlocarna option width in probabilistic mode
* Fix dependency in benchmark makefile
* Fix behavior of verbose mode / printing of consensus structure
* fix and support png/svg output in benchmarking
* Fix: use reference type for loop var
* Fix bug: quote names in newick trees
* Fix buffer size

### Internal changes:

* Update configure.ac by autoupdate
* Strip -fno-lto from Vienna RNA configured VRNA_LIBS
* doxygen: Turn off autolink support
* Require Vienna RNA 2.5.1
* Update tests

---

## 2.0.0RC9

### MLocarna functionality

* The calculation of the guide tree is rewritten; in consequence mlocarna is now more useful
    for clustering RNAs, since score distances are reported in the tree
    output (as well as names that correspond to intermediate alignments of
    subsets of the input sequences in subtrees)

* Generally, improved functionality for multiple alignment

    - better support for 'clustering' use case
        - distances and inner node labels in newick output
        - inner node labels correspond to names of intermediary alignments
    - introduces new perl Tree class with internal full tree representation
    - supports flexible traversals, filtering, etc
    - parsing and writing newick format fully supported
    - handles multi-ary trees
    - assign meaningful indices to intermediary alignments from most closely related to further apart
    - changes naming scheme of intermediary alignments in iterative
      alignment (include round in name)

* Perform distance-based upgma for guide tree computation:
    Calculate distances from locarna scores as
      distance  = max_score - score

* Add option (max-alignment-size) to limit the size of computed multiple alignments

* Don't write am and bm probabilities by default; introduce
options to control writing.

  Motivation: probabilities require lots of disk space, which is rarely
  used for inspection or further processing

  Introduces new options
    * --write-probs / --no-write-probs
    * --write-am-probs / --no-write-am-probs
    * --write-bm-probs / --no-write-bm-probs

  All of them default to 'no-write'

### Documentation improvements

* Improve man pages (add dois, slightly edit)
* Include disclaimer about the use of many threads in mlocarna man page
* Improve documentation of small tools (add missing man page includes)
* Add locarna_mcc man page

### Fixes

* Fix quoting in newick trees
* Update script mlocarna_nnames

### MLocarna minor or internal changes

* Refactorize guide tree computation.
* Inform user that results were written to target directory
* Code cleanup:
      Cleanup reference variables
      Cleanup amprobs/bmprobs hash references
      Cleanup sparse matrix module, return hash refs
      Change global var %opts to hash reference
      Expand leading tabs
      Cleanup in MatchProbs module
      Cleanup threaded/unthreaded modules
      Cleanup match probability computation.
      Remove redundant code for sequential computations of already parallelized steps

* Avoid unnecessary cloning of hashs; cleanup shared cloning
* factorize: introduce name normalizer class

### Other internal changes

* Remove obsolete script to label trees
* Update CARNA server related scripts
* Slightly simplify make of man pages (require man_include files)

## 2.0.0RC8

### Bugfixes

* Fix bug in sequence trace probability computation, which was
    potentially critical for large instances (e.g. aligning sequences with ~1000
    length difference)

* Remove superfluous output of alignments in probabilistic mode

### Code improvements/Maintenance

* Update travis CI script
* Fix initialization order [FitOnOff]
* Fix missing and superfluous headers [alphabet.hh/icc]

### Experimental Feature

* experimental extension of boundary prediction by position wise
    penalties, can use consistency file information [reliability-profile.pl, locarnap_fit]


---


## 2.0.0RC7

* Accept empty sequences in input
* Re-enable computing consensus structures in mlocarna output and stk-files
* Remove perl prototype specs to avoid perl 5.18 incompatibility

## 2.0.0RC6

* Adjust behavior of locarna and mlocarna for pairwise alignments
    - do not use alifold in binaries when calculating base pairs for single sequences
    - use noLP in mlocarna for calculating base pair probabilties /and/ the alignment
    - increase output precision for probabilities in pp files
       (reducing rounding issues for the typical cutoff)

## 2.0.0RC5 (2018-05-08)

* Allow penalized in global alignment

## 2.0.0RC4 (2018-05-04)

### Bugfixes

* Fix kbest HIT line
* Fix static linking (with gomp)
* Remove one extra line after CLUSTAL header (which used to be written by mlocarna)

## 2.0.0RC3 (2018-04-30)
### Bugfixes

* Fix kbest alignment (allow and detect empty alignments)
* Fix local end output (HIT line) for non-local alignments

### Features / Improvements

* Add stopwatch (for consensus computation)
* Introduce parameter --consensus-gamma for mea consensus 'gamma'

### API

* Class Alignment handles start and end position for non-local
    alignments; methods local_start(), local_end() of Alignment are
    replaced by start_positions(), end_positions()

## 2.0.0RC2 (2018-04-16)
### Bugfixes

*   Fix free endgaps bug
         (negative alignment scores were handled incorrectly)
*   Fix empty anchor constraint
         (correctly recognize/report when anchors are specified)

### Features / Improvements

Both improvements were introduced to handle large instances

* Improve alifold consensus computation
        - turn off by default (also speeds up mlocarna)
        - handle cases where consensus is not/cannot be computed
        - warn if consensus is not computed / fail if explicitely
          required (pp)
        - handle recomputation of consensus if file/screen output differ
          (local out)
* Support sequence envelope computation with extended precision
        - set the default to extended precision
        - recognize failing partition function calculation
        -   (for sequence trace and match probabilities)
        - extend LocARNA's quadmath.hh header and simplify its use
        - use templates for partition function type
        - support extended (long double) and quad precision (float128)


## 2.0.0RC1 (2017-12-13)
### Features / functionality

* Sequence envelopes for locarna, locarna_p and sparse;
    - Compatible to anchors, local alignment, and free endgaps
    - New parameter: --min-trace-probability (on by default)
    - Refactor class MatchProbs -> EdgeProbs and sub-classes;
      improve interface and docu
    - use sequence envelopes by default

* Change default parameters for scoring (new gap cost and tau factor)

* Support higher precision floating point arithmetic for locarna-p
    - pf calculation (long double and quad math) selectable dynamically without recompiling
    - support extended (long double) and quad precision (quadmath)

* Support caching of pp files (--dp-cache)

* Support in-package benchmarking (measuring sps, sci, and mcc)

### Security
* Avoid shell evaluation in mlocarna system calls

### Optimizations and code improvements

* Modernize and restructure code
        Transition to C++14 standard:
            use of unique_ptr, auto types,
            range loops, template meta programmin using type_traits, ...
        Systematic clean up of the library classes
        Introduce named arguments; zip and enumerate iterators
        Improved organization of command line help; reduce redundancies
        Improve handling of format and annotation types in MultipleAlignment

* Optimize inner loop: add sentinels to left and right adjacency lists
* Reduce memory footprint of locarna

* Anchor constraints restrict trace controller

* Refactor computation of input pp files [mlocarna]
* Use locarna_rnafold_pp for constrained folding [mlocarna]

* Update travis CI test

### Bugfixes

* Fix specifying anchor constraints via bed file [mlocarna]
     allow arbitrary order of entries in bed input
* Fix ribofit code (make experimental ribofit scoring usable)
* Fix confusion matrix: count_common_bps()

---


## 1.9.2 (2017-07-04)
    Improve consensus dotplot computation (avoid accumulation of low
      probability base pairs)
    Some cleanup

    1.9.2pre4
    Fix left and rightmost anchors (fixed local alignment); add unit tests
    Optimize maximization in innermost loop
    Add unit and program tests

    1.9.2pre3
    Make consensus structure / dotplot types independent [option
      --consensus-structure alifold (current default) used to also
      determine the computation of consensus dotplots. The consensus
      dot plot is now independent (new default: no alifold consensus
      dot plot).]

    fix strict semantics for anchor constraints
    fix compatibility '--keep-sequence-order' with anchor constraints

    1.9.2pre2
    mlocarna: write ps files (aln.ps, alirna.ps and if probabilistic,
      alidot.ps) from alifold call on result alignment to the
      results output directory;
      use alifold options --aln --color --mis
    fix mlocarna: find pw-aligner and pw-aligner-p executables in PATH

## 1.9.1 (2017-01-22)
    improve help and man pages for locarna, locarna_p, sparse;
      refactor code for help text (string map for help text)
      write text inserts for man pages, which are autogenerated from help
    fix bug in reading stockholm anchors
    remove parameters min_am_prob and min_bm_prob from locarna
      and sparse CLIs;
      they are useful only for locarna_p
    introduce strict semantics for anchor constraints
      and make default (see mlocarna --relaxed-anchors)
    fix exparna_p output
    add mlocarna option to realign; support 'consensus'
      constraint specifications


## 1.9.0 (2016-Dec-26)
    add mlocarna option to load (default) configuration from file
    add mlocarna option to read anchor constraints from
      bed-format specification
    add option --maxBPspan (analogous to ViennaRNA's maxBPspan)
    activate penalized option for mlocarna
    update documentation of mlocarna
    update reliability-profile.pl: suppress too verbose output;
      improve behavior of subseq option
    Tests:
        introduce mlocarna test with comparison to reference results
          - generate reference results with make gen-test-results
        add stockholm writing test
    Bug fixes:
        catch empty sequence case for alifold
        fix error message concerning max-diff-aln in mlocarna
        fix mlocarna bug: write stockholm also in probabilistic mode
        fix mlocarna error handling when running multithreaded
        fix doxygen/docu generation


## 1.8.12 (2016-Dec-1)
    Add checks: multi-threaded mlocarna, locarnate
    Fix severe problem with threads in mlocarna, which caused hangs
    Fix mlocarna with sparse (precomputation of inloop probs)
    Update VRNA/C11 configuration (adapts automatically; remove --enable-STDC )
    Require ViennaRNA >=2.3.1
    Improved options help
    Cleanup some files in the repository and distribution
## 1.8.11 (2016-Aug-15)
    Add configure flag --enable-STDC; set default to C11 (RNAlib compatibility)
    Fix output order of reliable structures in mlocarna
    Migrate to new vrna api; fix pf_scale bug
## 1.8.10  (2016-Jul-15)
    Improve input checking in mlocarna
    Increase API version to 1.5
    Wrap long alifold consensus output in mlocarna
    Add consensus structure in locarna output, selectable types mea|alifold|none
    Add optional stockholm output to mlocarna
    Introduce quiet switch for locarna, locarna_p, sparse, and exparna_p
    Code cleanup in driver programs
    Handle sequence name clashes by systematic renaming
    Read stockholm format
    Factorize common code in locarna, locarna_p and sparse (MainHelper)
    Add stockholm support input and output to locarna
    Fix bug occuring with empty input sequences
    Call alifold with /best/ parameters (from perl /and/ c++)
    Fix locarnate (compatibility to current locarna); add unit test

## 1.8.9   (2016-21-04)
    improve configuration
    improve ubuntu packaging
    fix anchor constraint check in mlocarna fasta input
    allow ribofit for mea - e.g. for probabilistic; this used to require ribosum

## 1.8.8   (2016-15-04)
    adapt to and require ViennaRNA >=2.2.4
    fix all static linking
    link binaries statically to libLocARNA by default;
    otherwise use --enable-dynamic-binaries
    fix & clean up code (after cppcheck static analysis)
    prepare for rpm packaging


---

# 2015

## 1.8.7    (2015-11-12)
    prepare for Ubuntu packaging
    fix man pages
      - update generation (add missing man pages)
      - fix program names
      - improve/fix help/man of binaries and locarnate
    clean exploc_p (get rid of "in"-mechanism)
    fix problems in tests (make check) w/o enable-debug
    remove wrong MultipleAlignment format check
    remove locarna.pl from distribution (functionally replaced by binary)
## 1.8.6   (2015-10-31)
    remove unused parameter "width" from locarna_rnafold_pp
    Increase API version to 1.4
    Get rid of locarna.hh
    In clustal(-like) alignment output, adapt formatting to long names
    Improve libRNA config, in cases where pkg-config does not work:
    allow compiling without libgomp, if libgomp is unavailable and
    is not used by libRNA
## 1.8.5a
    Fix behavior in case structure constraints + local folding (see 1.8.2)

## 1.8.5   (2015-10-10)
    allow empty input sequences/alignments and fix related problems
    fix score of local alignment w/ empty input sequences

## 1.8.4   (2015-10-6)
    fix potential overwriting of parameters in classes Aligner*,
    which was observed only with gcc-5.2
    fix bug that crashed printing string options in verbose mode
## 1.8.3
    improve portability for easier installation on recent MacOS X:
    compile without openmp; handle VRNA w/o pkgconf file
    fix c++11 warnings; fix incompatibility with libc++
    fix incompatibility with clang (now used by MacOS X); fix warnings of clang
    remove g++ warnings
## 1.8.2
    change semantics w.r.t. constraints + local folding:
        constraints turn off local folding
    reformat change log
    support generating doxygen documentation from
        distribution by make doxygen-doc
    penalized alignment
    experimental galaxy wrapper generation
## 1.8.1
    implement ribofit
    Bugfix in OMatrix: incorrect address calculation
    Bugfix: LocARNA crashed while writing the result pp file if
    alignments of zero length are computed with options local file
    output and alifold consensus dp (alifold workaround)
    bugfixes (OMatrix, alifold len 0 bug)
    implement slide rule and conflict rule for f1 and mcc calculation;
    also use default minimal base pair length 3 in class
    ConfusionMatrix; class interface changes slightly


---

# 2014

## 1.8.0
    Remove absolute path dependencies in scripts. This changes the
    behavior, since paths are not fixed at compile time. For
    example, the Vienna RNA programs are located via PATH environment
    variable.
    Require linking to VRNA. Remove configure options. Disallow
    renaming of Vienna RNA tools RNAfold, etc.; expect them in path.
    Always compile SPARSE/ExpaRNA-P. Don't pass mlocarna script
    through configure anymore.
    Add option --sparse to mlocarna
    Add rnafold-temperature option to mlocarna
    Use rnafold temperature and parameter files for alifold
    Warn if used with alifold-consensus-dp
    Add unit tests for mlocarna calls (standard, prob, sparse)
    Align alifold structure to alignment in mlocarna output
    Add Example gcvT.fa for SPARSE
    sparse: --special-gap-symbols option added (default is false)

## 1.7.17
    max-bps/uil/bpil-length-ratio for SPARSE/ExparnaP
    Second type of stacking terms: new-stacking
    Add heuristic to restrict maximum number of "significant" base
    pairs: max-bps-length-ratio
    Add max-diff restriction at arc match ends: --max-diff-at-am
    Restructure parameter passing to Aligner* classes: use "named
    argument" idiom
## 1.7.16
    ExpaRNA-P rename of tools
## 1.7.15
    Fix bug concerning iterative alignment and alignment evaluation
    Rename locarna_X to exparna_p
    ExpaRNA-P code "publication ready"

## 1.7.14
    Handle inconsistencies due to anchor constraints in conflict
    with other constraints/heuristics, e.g. max-diff: 1) accept -inf
    pairwise similarities for guide tree construction; 2) clear
    consensus anchor strings with duplicate names



---

# 2013

## Jan 2013
    Fix bug. Quotes in fasta-names caused failure when
    generating pps via locarna_rnafold_pp in mlocarna.
    (1.7.3.1)

    Resolve rnalib dependency. New configure option with-sparse
    is now required to compile with ExpARNA-P and SPARSE functionality.
    (1.7.5)

## Mar 2013
    Merge in new sparsification related code.
    Added fasta input check in mlocarna.
    (1.7.6)
    Small bug fixes
    Fix memory leaks
    Prepare alignment evaluation functionality
    Code cleanup/refactorizations
    (1.7.6.9)
    Introduce configuration header config.h
    Fix man page dependency problem

## May/Jun 2013
    Make sure that make check and make distcheck work
    Split RnaData into RnaData and RnaEnsemble
    Introduce new version of pp format (2.0)
    Read/write of inloop probabilities
    Cleaner separation of reading/writing of alignments/probability data v.s. computation;
      in particulare remove I/O from Alignment class
    Cleaner handling of gaps, alignment edges, and gap symbols in I/O
    New option --stopwatch (in locarna*) to turn on detailed internal timing
    Change class Sequence to be a mere interface to MultipleAlignment that allows operator [] for column access
    Make constraints optional part of multiple alignments; make write/readable with the alignment in clustal and pp format
    Remove options anchorA, anchorB from locarna and locarna_X
    Support anchor and structure annotation in MultipleAlignment
    RnaEnsemble constrains the folding by structure annotation of
    Multiple Alignment (if given)
    (Release 1.7.7)
## Jul 2013
    Support precomputing in loop probabilities for SPARSE
    Fixed reading of PP 2.0
    (Release 1.7.8)
## Sep 2013
    LocARNATE: new code, improvements
    Exparna-P: Suboptimal traceback and improvements
    Change default behavior of mlocarna for sequence-local alignment: proceed
    with only-local subalignments; options local-progressive and
    global-progressive control behavior, where --global-progressive
    reactivates the previous one.
    (Release 1.7.9)
## Oct 2013
    Newick bugfix
    Small modifications for library use (e.g. for carna)
    (Release 1.7.10)
## Nov 2013
    Improve max-diff heuristic (without reference alignment): make
    symmetric and always allow sufficiently large deviation
    (Release 1.7.11)
    Fix bug in locarna_rnafold_pp (reading clustalw from stdin)
    Fix bug concerning structure constraints
    Update format conversions in mlocarna to pp 2.0 format
    Support for fixed structure constraints (directly in C++)
    (Release 1.7.12)
    Fix bug that compromised LocARNA-P w/ anchor-constraints in 1.7.12
    In-loop support for fixed structure constraints
    (Release 1.7.13)
    Bugfix: computation of consensus anchor constraints could fail; now resolve
    name clashs by selecting lexicographically smaller name.
    (Release 1.7.13.1)


---

# 2012
## Jan 2012
    fix bugs in reliability plotting and "fitting"
    (1.6.2.3)
    New tool locarna-mea for computing the max. expected accuracy structure
    from a dot plot
    (1.6.2.4)
    New features for pp2dot script
    New script for generating average pp-files
    In combination, generate nice colored dot plots;
    used first for the output of the carna server
    (1.6.3)
## Mar 2012
    Merge with recent Exparna-P/locarna_X code
    (1.7)
    New classes for infinite integers (this is used for scores,
    which can be -infty, i.e. undefined). The new mechanism automatizes
    normalization of infinity using the type system and avoids wrong
    usage by the type system. Consequently, it is more flexible. It avoids
    a bug in normalized alignment and helps to clean up the Exparna-P code.
    (1.7.1)
## Apr/May 2012
    Add reading of score matrix in mlocarna
    Fix bug in probabilistic mode, which did not allow '#' in sequence names
    Fix bugs in output of score matrix
    Fix drawdplin
## June 2012
    Fix bug in trace controller
    (1.7.2)
    Add locarna_n code, alpha status
## Aug 2012
    Fix a bug related to use of trace controller (max-diff) with free
    endgaps.
    Fix autodetection of input multiple alignment file format,
    i.e. FASTA vs. CLUSTAL.
    Change interface of scoring such that gapA (gapB) depends on only
    positions in sequence A (B), resp.
    (1.7.2.1)
    Major improvements of the RnaData class and its interface
    Some clean up
    (1.7.2.5)
    Some bug fixes, support for graph kernel guide tree
    (1.7.2.6)
    Include graph kernel guide tree support in mlocarna
    (1.7.2.7)
    Some cleanup, release.
    Fix bug (missing initialization) in RnaData
    (1.7.3)



---

# 2011

## Jan 2011
    Alignment in limited deviation to a reference alignment
    program deviation computes deviation of alignment to reference
    Deactivate max-diff-aln for --probabilistic/locarna_p
    (Version 1.6)

## Feb 2011
    Restructure code
    library libLocARNA is generated and is used as library by executables

    Fix name normalization bug in mlocarna.in
    (Version 1.6.1)
## Mar 2011
    Doxygen compiles without warning and produces useful documentation
    of the library

    Vienna RNA library is now required for compilation and will be linked

    Locarna accepts clustal files as input and then computes base pair
    probabilities via pf_fold of the Vienna RNA lib. Currently this is
    limited to single sequence clustal files.
    Implement reading of mfasta files (still switched off).

## Apr 2011
    Implement several alignment comparison measures (and use in locarna_deviation) (EXPERIMENTAL)

    Make linking to Vienna librna configurable, default=off! Activate
    by calling ./configure with --enable-librna.

    New strategy for consensus dot plot generation: employ RNAalifold
    -p.  Enable via mlocarna's option --alifold-consensus-dp Note: the
    default consensus dot plot computation seems to be flawed! The new
    tool pp2dot allows to visualize consensus dot plots (by converting
    pp to Vienna-like dot plot postscript files) and should help to
    improve the default strategy (otherwise we are going to replace
    it!)

    mlocarna accepts via --tree-file a super-tree as guide tree and
    automatically projects to the taxa in the input.

    Allow quotation in newick trees (using Text::ParseWords)
    and use quotations when generating newick trees

## May 2011
    Fix a problem with compiling on case-insensitive file systems (generate locarna.bin instead of locarna executable in src, rename at install)

    add alifold-consensus-dp option for mlocarna

## Aug 2011: Add code for 'locarna_X'

## Sep 2011: Adding Exparna-P code for "in loop" probabilities

## Oct 2011: pre-release 1.6.2

## Nov 2011:
    Add mlocarna options (pw-aligner and pw-aligner-p) to run
    specified binaries in place of locarna and locarna_p. This can be
    used to plug in compatible programs like carna

    Bug fix; compilation without linking libRNA was broken.
    release (1.6.2)

## Dec 2011
    Bug-fixes: treefile regexp, LARGE_PF name clash mit Vienna package
    (1.6.2.1)
    add script to translate to normalized names for mlocarna

    Add options pw-aligner-options and pw-aligner-p-options, to mlocarna to pass option strings to the pairwise
    aligner.  This is useful in combination with specifying
    non-standard pairwise alignment programs.

    Support generation of pp files from dp2.ps and dp.ps files when
    skipping dot plot computation with --skip-pp (keep name for
    legacy).  This allows the combine user-defined dotplots with
    anchor constraints from the mfasta input file.

    Support fixed structures (with pseudoknots) in the fasta input
    of mlocarna as structure constraint string with tag #FS
    All base pairs in the fix structure receive a probability of 1.
    (1.6.2.2)


---

# 2010

## Jan 2010
    new tool Utils/reliability-profile.pl generates nice
    profile plots from mlocarna --probabilistic output directories

    stacking turned off, unless exp-prob option used

    (locarna-1.5a)

## Mar 2010
    Add normalized local alignment, proof of concept
    (locarna-1.5.1)
    Combine normalized alignment and kbest
    Introduce option better to enumerate suboptimal solution
    better than given threshold
    --kbest and --better can be used together
    (locarna-1.5.2)
    Fix build in separate build directory

## May 2010
    Introduce new representation of fasta data
    Major restructuring/refactorization of mlocarna.in and MLocarna.pm
    Move many functions of mlocarna.in to packages and
    introduce several new packages to group functions thematically
    (locarna-1.5.3)

    Resolve a bug that made LocARNA handle Ts as unknown
    characters. LocARNA now automatically transforms Ts to Us in input
    sequences. The bug did not affect input sequences without Ts. For
    matchs involving Ts, LocARNA wrongly used an ad-hoc score instead
    of the ribosum score for matches involving 'U', as one would
    expect.

## Oct 2010
    Add EdgeController: adds feature to constrain (max-diff-match)
    alignment relative to a given alignment.  Change handling of
    max-diff-match in case of unequal length: instead of relaxing
    constraint, change constraint to |i*(lenB/lenA)-j|<=max-diff
    Add multiple alignment class to support construction of EdgeController.

## Nov 2010
    Change Copyright to GPL v2 and update README

    rename max-diff-match to max-diff heuristic.

    replace EdgeController by TraceController use the max-diff delta
    ( TraceController ) to restrict alignment of paired and unpaired
    positions (not only ends of arcs) in Aligner. Semantics is given by
    restriction of alignment cuts/matrix cells!


---

# 2009

## Jan 2009 (v1.42)
    iterative alignment using SCI as optimization criterion

## Feb 2009 (v1.4.4alpha)
    cleanup mlocarna code
    cleanup mlocarna output (inclusive output modes -v, --moreverbose)
    improve reliability computation
    multi-threaded mlocarna, new option --cpus <n>: use n threads
    more structure for mlocarna output to file system (distribute
    files to subdirs)
    introduce --pf-scale parameter to locarna_p and mlocarna

    precompute exp base match similarity and gap cost for partition
    function => significant speed up of locarna_p
    change handling of max_diff: guarantee at least |lenA-lenB|
    improve outside algorithm in locarna_p (reduce sizes of Mrev-matrices)
    => with prior optimization massive speedup of locarna_p

    add probability filtering after cbt
    further improvement of outside (cheaper computation of certain
    Mprime entries)
    iteration with reliable structure as constraints
    (option --it-reliable-structure)

    support for non-threaded perl

## Mar 2009
    fix some problems with the option parsing
    (e.g. now -p=0.1 is correct and -p=a throws an error)
    fix bug in initialization for structure locality
    change semantics of constraints in the context of locality:
    anchor constraints (i,j) now enforce that edge (i,j) is part
    of the local alignment (for sequence and structure locality)
## Apr 2009
    add missing options to mlocarna, fix bug in locarnate (thanks to Wolfgang)
    fix bug with sharing/threads in mlocarna (v1.4.6)
    locarna_p computes probabilities for matching fragments
## May 2009
    clean manpage of mlocarna
## Jun 2009
    offer large data type "long double" for partition function
    computation.  this is available when configuring through option
    --enable-large-pf. Use this in case of long alignments, where the
    partition function exceeds the range of type "double". Although
    --pf-scale sometimes provides a workaround, for multiple alignment
    of large sequences it is often the only viable choice.

## Jul 2009
    add free endgaps via option --free-endgaps "----", for locarna
    gaps at all for ends can be allowed for free independently.
    For mlocarna, --free-endgaps will allow free end gaps at any end.

## Sep 2009
    fixed serious bugs concerning stacking: due to these bugs,
    mlocarna did not profit from stacking and locarna could get wrong
    results with stacking turned on.

    Change base pair filtering by probability for stacking: ignore
    stacking probability, filter only by base pair probability.

    (locarna-1.4.9alpha)


---

# 2008
## Feb 2008
    remove bug in anchor constraints
    Ribosum-85-60 compiled in
    changed default parameters to something reasonable

## Mar 2008
    remove bug in rna_data.cc reading dp.ps, v 1.2.4

## Jun 2008
    change char * to std::string in options module
    (the g++-compiler does not like it anymore :))

## Jul 2008
    integrate locarnate - a tool for multiple alignment
    that uses locarna in combination with tcoffee
    (Version 1.3)

## Aug 2008
    partition function version of locarna (locarna_p)
    experimental support for multiple alignment using
      basematch probabilities and probabilistic
      consistency transformation
    (Version 1.3.2)

    mlocarna can write basematch probabilities for plotting
    plotting of basematch probs (Utils/plot-bmprobs)

    mlocarna now correctly accepts input trees in newick-format that
    are terminated by ";". However, it does not enforce the ";"-termination
    mlocarna now writes the result.tree in ";"-terminated newick-format

    use rounding (instead of truncate) in scoring module, when
    converting from double to score_t. This fixes some inconsistency
    between 32/64bit and is the correct thing to do anyway

## Oct 2008
    improve speed up gained from using anchor constraints

## Nov 2008 (v1.4)
    add probabilistic scoring with consistency transformation
    with structural scoring in progressive alignment steps

    change options of mlocarna for probabilistic/mea-type alignment



---

# 2007
## Sept 2007 (1.0alpha)
    refactoring; new classes Scoring, RnaData, ...
    preparations for better scoring:
    reading Ribosum-Data, new class Ribosum
    computing match probabilities in sequence alignment (Gotoh):
        Implemented methods: pairHMM, partition function, local pf
        new class MatchProbs
    improve option parser: introduce sections and section hiding
    integrate new scoring scheme: maximum expected accuracy
    significant speed-up due to re-organization of the recursion
    equations:
       redefine D to already contain the score contribution of
       the arcmatch. This reduces these computations by orders of
       magnitude!

## Nov 2007
    add constraints to pairwise and multiple alignment. This means
    for locarna: integrate anchor constraints
    for mlocarna: read and pass structure constraints to RNAfold -p,
                  read and pass sequence constraints to locarna
    The pp-format is redefined in order to pass constraint information from
    mlocarna to locarna.
    New class for representing and calculating with infinite scores
    changed default sequ-local=on to off
    increased min_prob in mlocarna.pl

    bug fix in computation of consensus probs

    rewrite noLP
    add (pseudo)affine gap penalties
    make mea working
    change default output width to 120
    ribosum scoring for classic score
    change handling of options in mlocarna

    mlocarna learns iterative refinement
    mlocarna learns conistency based transformation

## Dec 2007
    true affine gap cost (implemented using additional vectors E and F)
    implemented stacking (v_1_2alpha)
    removed bugs in traceback
    cleaned options
    include version name in output (v_1_2_1a)
    add man pages
    add Ribosum85-60 matrix in distribution
    avoid writing to global /tmp directory (v_1_2_2a)
    avoid potential file name conflicts in mlocarna
    accept duplicate (or not prefix unique) names in mfasta input


---

# 2006

## Jan 2006:
    add support for aligning alignments
    (used for multiple progressive alignment)

## Aug 2006:
    Lower memory consumption of locarna.

    In particular important in case of
    scanning (target against genome)

    Add some options/features that are useful for scanning
    option --kbest
    option --max-diff-am
    option --ungappedA
    option --local-output
    handling of nucleotide-symbol N for unknown in input


## Oct 2006
    add support for stacking alignment, not finished
## May 2006 (0.991)
    small changes, rename mlocarna-p to mlocarna


---

# 2005

## Mai 2005
    global alignment
## June 2005
    extension to local alignment (sequence+structure locality),
    no_lonely_pairs option
## Nov 2005
    command line options
## Dez 2005
    recognize and read dot.ps files (in add to formerly used pp-files)
