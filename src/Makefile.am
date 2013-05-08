## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


########################################
### Rules for LocARNA library libLocARNA

lib_LTLIBRARIES=libLocARNA-@API_VERSION@.la

libLocARNA_@API_VERSION@_la_SOURCES=LocARNA/aux.cc			\
	LocARNA/sequence.cc LocARNA/stopwatch.cc LocARNA/options.cc	\
	LocARNA/basepairs.cc LocARNA/alignment.cc LocARNA/rna_data.cc	\
	LocARNA/scoring.cc LocARNA/infty_int.cc				\
	LocARNA/match_probs.cc LocARNA/ribosum.cc			\
	LocARNA/anchor_constraints.cc LocARNA/arc_matches.cc		\
	LocARNA/aligner.cc LocARNA/stral_score.cc			\
	LocARNA/aligner_p.cc LocARNA/trace_controller.cc		\
	LocARNA/multiple_alignment.cc LocARNA/rna_structure.cc		\
	LocARNA/confusion_matrix.cc

libLocARNA_@API_VERSION@_la_LDFLAGS = -version-info $(SO_VERSION)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ../LocARNA-$(API_VERSION).pc

library_includedir=$(includedir)/LocARNA-$(API_VERSION)
nobase_library_include_HEADERS = locarna.hh LocARNA/aux.hh		\
	LocARNA/stopwatch.hh LocARNA/options.hh LocARNA/matrices.hh	\
	LocARNA/sparse_matrix.hh LocARNA/sequence.hh			\
	LocARNA/basepairs.hh LocARNA/alignment.hh LocARNA/rna_data.hh	\
	LocARNA/params.hh LocARNA/scoring.hh LocARNA/infty_int.hh	\
	LocARNA/match_probs.hh LocARNA/alphabet.hh			\
	LocARNA/alphabet.icc LocARNA/ribosum.hh				\
	LocARNA/anchor_constraints.hh LocARNA/arc_matches.hh		\
	LocARNA/aligner_restriction.hh LocARNA/aligner.hh		\
	LocARNA/aligner_p.hh LocARNA/stral_score.hh			\
	LocARNA/trace_controller.hh LocARNA/multiple_alignment.hh	\
	LocARNA/ribosum85_60.icc LocARNA/rna_structure.hh		\
	LocARNA/confusion_matrix.hh LocARNA/discrete_distribution.hh	\
	LocARNA/string1.hh LocARNA/type_wrapper.hh

## binary programs
##
##   use extension .bin for binary locarna to avoid name collission in src dir
##
bin_PROGRAMS = locarna.bin ribosum2cc locarna_p locarnap_fit	\
               locarna_deviation locarna_rnafold_pp


if HAVE_LIBRNA
    libLocARNA_@API_VERSION@_la_SOURCES += \
        LocARNA/mcc_matrices.cc
    nobase_library_include_HEADERS += \
        LocARNA/mcc_matrices.hh
endif


## whether to build components that depend on strong 
## probabilistic sparsification via the sparsification mapper
## requiring direct access to the McCaskill matrices
if BUILD_SPARSE
    libLocARNA_@API_VERSION@_la_SOURCES += \
	LocARNA/aligner_n.cc LocARNA/sparsification_mapper.cc \
	LocARNA/exact_matcher.cc

    nobase_library_include_HEADERS += \
	LocARNA/aligner_n.hh LocARNA/sparsification_mapper.hh \
        LocARNA/exact_matcher.hh

    bin_PROGRAMS += locarna_X locarna_n
endif

#remove the extension .bin for installation
transform = s/\.bin$$//

# all binaries should be linked to libLocARNA
LDADD=libLocARNA-$(API_VERSION).la

BUILT_SOURCES=LocARNA/ribosum85_60.icc

check_PROGRAMS = test_trace_controller test_rna_data
test_trace_controller_SOURCES = test_trace_controller.cc
test_rna_data_SOURCES = test_rna_data.cc
TESTS = $(check_PROGRAMS)

dist_locarna_bin_SOURCES = locarna.cc

locarna_p_SOURCES = locarna_p.cc

locarna_X_SOURCES = locarna_X.cc 

locarna_n_SOURCES = locarna_n.cc


nodist_locarna_bin_SOURCES = LocARNA/ribosum85_60.icc

ribosum2cc_SOURCES = ribosum2cc.cc

locarnap_fit_SOURCES = fit.cc LocARNA/fitonoff.cc LocARNA/options.hh LocARNA/fitonoff.hh

locarna_deviation_SOURCES = deviation.cc

locarna_rnafold_pp_SOURCES = locarna_rnafold_pp.cc


LocARNA/ribosum85_60.icc: ribosum2cc
	./ribosum2cc Ribosum85_60 "$(top_srcdir)/Matrices/RIBOSUM85_60" >$@

CLEANFILES=LocARNA/ribosum85_60.icc
