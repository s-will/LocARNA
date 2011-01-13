bin_PROGRAMS = locarna ribosum2cc locarna_p fit locarna_X deviation

BUILT_SOURCES=ribosum85_60.icc


MYSOURCES=options.cc options.hh \
	aux.hh aux.cc \
	sequence.cc sequence.hh basepairs.cc basepairs.hh \
	alignment.hh alignment.cc aligner.hh aligner.cc \
	params.hh matrices.hh rna_data.cc rna_data.hh \
	sparse_matrix.hh scoring.hh scoring.cc \
	match_probs.cc match_probs.hh alphabet.hh alphabet.icc \
	ribosum.hh ribosum.cc anchor_constraints.hh anchor_constraints.cc \
	infty_arith_int.hh infty_arith_int.cc \
	stral_score.hh stral_score.cc \
	arc_matches.cc arc_matches.hh \
	aligner_p.cc aligner_p.hh \
	trace_controller.cc trace_controller.hh \
	evaluator.cc evaluator.hh \
	multiple_alignment.hh multiple_alignment.cc

check_PROGRAMS = trace_controller_test
trace_controller_test_SOURCES = trace_controller_test.cc $(MYSOURCES)
TESTS = $(check_PROGRAMS)


dist_locarna_SOURCES = locarna.cc $(MYSOURCES)

locarna_p_SOURCES = locarna_p.cc $(MYSOURCES)

locarna_X_SOURCES = locarna_X.cc $(MYSOURCES) \
	exact_matcher.hh exact_matcher.cc \
	basepairs_looptraversal.hh basepairs_looptraversal.cc

nodist_locarna_SOURCES = ribosum85_60.icc

ribosum2cc_SOURCES = ribosum2cc.cc ribosum.cc ribosum.hh matrices.hh alphabet.icc alphabet.hh

fit_SOURCES = fit.cc fitonoff.cc fitonoff.hh options.hh options.cc

deviation_SOURCES = deviation.cc multiple_alignment.hh multiple_alignment.cc aux.hh aux.cc

ribosum85_60.icc: ribosum2cc
	./ribosum2cc Ribosum85_60 $(top_srcdir)/Matrices/RIBOSUM85_60 >$@

CLEANFILES=ribosum85_60.icc
