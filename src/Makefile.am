SUBDIRS=Utils . Tests

AM_DEFAULT_SOURCE_EXT = .cc

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
	LocARNA/basepairs.cc LocARNA/alignment.cc			\
	LocARNA/rna_ensemble.cc LocARNA/rna_data.cc			\
	LocARNA/pfold_params.cc LocARNA/scoring.cc			\
	LocARNA/infty_int.cc LocARNA/match_probs.cc			\
	LocARNA/ribosum.cc LocARNA/ribofit.cc				\
	LocARNA/sequence_annotation.cc LocARNA/anchor_constraints.cc	\
	LocARNA/arc_matches.cc LocARNA/aligner.cc			\
	LocARNA/stral_score.cc LocARNA/aligner_p.cc			\
	LocARNA/trace_controller.cc LocARNA/multiple_alignment.cc	\
	LocARNA/rna_structure.cc LocARNA/confusion_matrix.cc		\
	LocARNA/global_stopwatch.cc LocARNA/mcc_matrices.cc		\
	LocARNA/aligner_n.cc LocARNA/sparsification_mapper.cc		\
	LocARNA/exact_matcher.cc LocARNA/params.cc

libLocARNA_@API_VERSION@_la_LDFLAGS = -version-info $(SO_VERSION)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ../LocARNA-$(API_VERSION).pc

library_includedir=$(includedir)/LocARNA-$(API_VERSION)
nobase_library_include_HEADERS = LocARNA/aux.hh LocARNA/plusvector.hh	\
	LocARNA/stopwatch.hh LocARNA/options.hh LocARNA/matrix.hh	\
	LocARNA/matrices.hh LocARNA/sparse_matrix.hh			\
	LocARNA/sparse_vector.hh LocARNA/sequence.hh			\
	LocARNA/basepairs.hh LocARNA/alignment.hh			\
	LocARNA/alignment_impl.hh LocARNA/rna_ensemble.hh		\
	LocARNA/rna_ensemble_impl.hh LocARNA/rna_data.hh		\
	LocARNA/ext_rna_data.hh LocARNA/rna_data_impl.hh		\
	LocARNA/ext_rna_data_impl.hh LocARNA/params.hh			\
	LocARNA/pfold_params.hh LocARNA/scoring.hh			\
	LocARNA/scoring_fwd.hh LocARNA/infty_int.hh			\
	LocARNA/match_probs.hh LocARNA/alphabet.hh			\
	LocARNA/alphabet.icc LocARNA/ribosum.hh LocARNA/ribofit.hh	\
	LocARNA/ribofit_will2014.ihh LocARNA/ribofit_will2014.icc	\
	LocARNA/sequence_annotation.hh LocARNA/anchor_constraints.hh	\
	LocARNA/arc_matches.hh LocARNA/aligner_restriction.hh		\
	LocARNA/aligner.hh LocARNA/aligner_impl.hh			\
	LocARNA/aligner_p.hh LocARNA/stral_score.hh			\
	LocARNA/trace_controller.hh LocARNA/multiple_alignment.hh	\
	LocARNA/rna_structure.hh LocARNA/confusion_matrix.hh		\
	LocARNA/discrete_distribution.hh LocARNA/string1.hh		\
	LocARNA/type_wrapper.hh LocARNA/global_stopwatch.hh		\
	LocARNA/tuples.hh LocARNA/mcc_matrices.hh			\
	LocARNA/aligner_n.hh LocARNA/sparsification_mapper.hh		\
	LocARNA/exact_matcher.hh LocARNA/main_helper.icc		\
	LocARNA/ribosum85_60.icc

## binary programs
##
##   use extension .bin for binary locarna to avoid name collission in src dir
##
bin_PROGRAMS = locarna.bin locarna_p locarnap_fit locarna_deviation	\
               locarna_rnafold_pp exparna_p sparse ribosum2cc

if STATIC_LIBLOCARNA
## link libLocARNA statically to the binaries
locarna_deviation_LDFLAGS=-static
exparna_p_LDFLAGS=-static
locarnap_fit_LDFLAGS=-static
locarna_bin_LDFLAGS=-static
locarna_p_LDFLAGS=-static
locarna_rnafold_pp_LDFLAGS=-static
ribosum2cc_LDFLAGS=-static
sparse_LDFLAGS=-static
endif

#remove the extension .bin for installation
transform = s/\.bin$$//

# all binaries should be linked to libLocARNA
LDADD=libLocARNA-$(API_VERSION).la

dist_locarna_bin_SOURCES = locarna.cc

locarna_p_SOURCES = locarna_p.cc

exparna_p_SOURCES = exparna_p.cc 

sparse_SOURCES = sparse.cc

locarnap_fit_SOURCES = fit.cc LocARNA/fitonoff.cc LocARNA/options.hh LocARNA/fitonoff.hh

locarna_deviation_SOURCES = deviation.cc

locarna_rnafold_pp_SOURCES = locarna_rnafold_pp.cc

BUILT_SOURCES = LocARNA/ribosum85_60.icc
CLEANFILES = LocARNA/ribosum85_60.icc

if !VC_CROSS_COMPILING
LocARNA/ribosum85_60.icc: \
	ribosum2cc \
	$(top_srcdir)/Data/Matrices/RIBOSUM85_60

	./ribosum2cc Ribosum85_60 "$(top_srcdir)/Data/Matrices/RIBOSUM85_60" $@
endif
