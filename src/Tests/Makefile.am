# all binaries should be linked to libLocARNA
LDADD=../libLocARNA-$(API_VERSION).la

AM_DEFAULT_SOURCE_EXT = .cc

BINTESTS = test_locarna_lib
SCRIPTTESTS = test_programs

test_locarna_lib_SOURCES = alphabet.cc anchor_constraints.cc		\
	catch.hpp ext_rna_data.cc matrices.cc multiple_alignment.cc	\
	rna_data.cc rna_ensemble.cc rna_structure.cc			\
	test_locarna_lib.cc trace_controller.cc zip.cc

TESTS= $(BINTESTS) $(SCRIPTTESTS)

check_PROGRAMS = $(BINTESTS)

MYTESTDATA    = archaea.aln archaea-aln.fa archaea-wrong.fa

MYTESTRESULTS = locarnate.testresult					\
	mlocarna-archaea-alifold.testresult				\
	mlocarna-bed-anchors.testresult					\
	mlocarna-haca-constraints.testresult				\
	mlocarna-probabilistic.testresult mlocarna-sparse.testresult	\
	mlocarna-threads.testresult mlocarna-realign.testresult		\
	exparna_p.testresult locarna-local.testresult			\
	locarna-normalized.testresult

BUILT_SOURCES = $(MYTESTDATA)

EXTRA_DIST    = $(MYTESTDATA:%=$(top_srcdir)/Data/Examples/%) \
	$(MYTESTRESULTS) \
	$(SCRIPTTESTS)

archaea.aln: $(top_srcdir)/Data/Examples/archaea.aln
	cp $< $@
archaea-aln.fa: $(top_srcdir)/Data/Examples/archaea-aln.fa
	cp $< $@
archaea-wrong.fa: $(top_srcdir)/Data/Examples/archaea-wrong.fa
	cp $< $@


CLEANFILES = $(MYTESTDATA)

## generate test results for mlocarna test
gen-test-results:
	srcdir=$(srcdir) \
	$(srcdir)/test_programs --generate-references
