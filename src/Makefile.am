AUTOMAKE_OPTIONS = subdir-objects

### Rules for LocARNA library libLocARNA

lib_LTLIBRARIES=libLocARNA.la
libLocARNA_la_SOURCES=LocARNA/aux.cc LocARNA/options.cc			\
	LocARNA/sequence.cc LocARNA/basepairs.cc LocARNA/alignment.cc	\
	LocARNA/rna_data.cc LocARNA/scoring.cc				\
	LocARNA/infty_arith_int.cc LocARNA/match_probs.cc		\
	LocARNA/ribosum.cc LocARNA/anchor_constraints.cc		\
	LocARNA/arc_matches.cc LocARNA/aligner.cc			\
	LocARNA/stral_score.cc LocARNA/aligner_p.cc			\
	LocARNA/aligner_p.hh LocARNA/trace_controller.cc		\
	LocARNA/evaluator.cc LocARNA/multiple_alignment.cc


library_includedir=$(includedir)/LocARNA
library_include_HEADERS = LocARNA/aux.hh LocARNA/options.hh		\
	LocARNA/matrices.hh LocARNA/sparse_matrix.hh			\
	LocARNA/sequence.hh LocARNA/basepairs.hh LocARNA/alignment.hh	\
	LocARNA/rna_data.hh LocARNA/params.hh LocARNA/scoring.hh	\
	LocARNA/infty_arith_int.hh LocARNA/match_probs.hh		\
	LocARNA/alphabet.hh LocARNA/alphabet.icc LocARNA/ribosum.hh	\
	LocARNA/anchor_constraints.hh LocARNA/arc_matches.hh
	LocARNA/aligner.hh LocARNA/stral_score.hh LocARNA/trace_controller.hh
	LocARNA/evaluator.hh LocARNA/multiple_alignment.hh

## binary programs
bin_PROGRAMS = locarna ribosum2cc locarna_p locarnap_fit locarna_deviation
# locarna_X 

# all binaries should be linked to libLocARNA
LDADD=libLocARNA.la

BUILT_SOURCES=ribosum85_60.icc

check_PROGRAMS = trace_controller_test
trace_controller_test_SOURCES = trace_controller_test.cc
TESTS = $(check_PROGRAMS)

dist_locarna_SOURCES = locarna.cc

locarna_p_SOURCES = locarna_p.cc

#locarna_X_SOURCES = locarna_X.cc \
#	exact_matcher.hh exact_matcher.cc \
#	basepairs_looptraversal.hh basepairs_looptraversal.cc

nodist_locarna_SOURCES = ribosum85_60.icc

ribosum2cc_SOURCES = ribosum2cc.cc

locarnap_fit_SOURCES = fit.cc fitonoff.cc LocARNA/options.hh fitonoff.hh

locarna_deviation_SOURCES = deviation.cc

ribosum85_60.icc: ribosum2cc
	./ribosum2cc Ribosum85_60 $(top_srcdir)/Matrices/RIBOSUM85_60 >$@

CLEANFILES=ribosum85_60.icc
