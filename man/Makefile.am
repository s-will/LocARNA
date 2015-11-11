HELP2MANOPTS=-Nl


dist_man_MANS=locarna.1 locarna_p.1 mlocarna.1 locarnate.1		\
         reliability-profile.pl.1 locarnap-predict-and-plot.pl.1	\
         locarnap-realign-all.pl.1 locarnap-revisit-RNAz-hits.pl.1	\
         locarnap-revcomp.pl.1 locarnap_fit.1 locarna_deviation.1	\
         exparna_p.1 exploc_p.1 sparse.1 locarna_rnafold_pp.1		\
         mlocarna_nnames.1

# dist_man_MANS=aln2fa.pl.1 alnsel.pl.1 aln-seqs.pl.1 average-dot.pl.1	\
# 	dot2pp.1 exparna_p.1 exploc_p.1 gen-reliab-dot.pl.1 locarna.1	\
# 	locarna_deviation.1 locarna-mea.pl.1 locarna-motif-scan.1	\
# 	locarna_p.1 locarnap_fit.1 locarna.pl.1				\
# 	locarnap-predict-and-plot.pl.1 locarnap-realign-all.pl.1	\
# 	locarnap-revcomp.pl.1 locarnap-revisit-RNAz-hits.pl.1		\
# 	locarna_rnafold_pp.1 locarnate.1 mlocarna.1	\
# 	mlocarna_nnames.1 plot-bmprobs.1 pp2dot.1			\
# 	reliability-profile.pl.1 ribosum2cc.1 sparse.1


## the man pages should *not* be remade when building unchanged distributions
## but should be rebuild after change of the options or package version.
## Thus, we make the depend on the sources and ensure that the binary
## is made before we call it (order prerequisite!)

# aln2fa.pl.1
# alnsel.pl.1
# aln-seqs.pl.1
# average-dot.pl.1
# dot2pp.1
# exparna_p.1
exparna_p.1: $(top_srcdir)/src/exparna_p.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/exparna_p
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/exparna_p

# exploc_p.1
exploc_p.1: $(top_srcdir)/Utils/exploc_p $(top_srcdir)/configure.ac
	-$(POD2MAN)  $(top_srcdir)/Utils/exploc_p >$@

# gen-reliab-dot.pl.1

# locarna.1
locarna.1: $(top_srcdir)/src/locarna.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarna.bin
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna.bin
#	strip .bin (special name handling for locarna binary)
	sed -i 's/locarna.bin/locarna/g' $@

# locarna_deviation.1
locarna_deviation.1: $(top_srcdir)/src/deviation.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarna_deviation
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna_deviation

# locarna-mea.pl.1
# locarna-motif-scan.1
# locarna_p.1
locarna_p.1: $(top_srcdir)/src/locarna_p.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarna_p
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna_p

# locarnap_fit.1
locarnap_fit.1: $(top_srcdir)/src/fit.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarnap_fit
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarnap_fit

# locarna.pl.1
# locarnap-predict-and-plot.pl.1
locarnap-predict-and-plot.pl.1: $(top_srcdir)/Utils/locarnap-predict-and-plot.pl \
	$(top_srcdir)/configure.ac
	-$(POD2MAN) $< >$@

# locarnap-realign-all.pl.1
locarnap-realign-all.pl.1: $(top_srcdir)/Utils/locarnap-realign-all.pl \
	$(top_srcdir)/configure.ac
	-$(POD2MAN) $< >$@

# locarnap-revcomp.pl.1
locarnap-revcomp.pl.1: $(top_srcdir)/Utils/locarnap-revcomp.pl $(top_srcdir)/configure.ac
	-$(POD2MAN) $< >$@

# locarnap-revisit-RNAz-hits.pl.1
locarnap-revisit-RNAz-hits.pl.1: $(top_srcdir)/Utils/locarnap-revisit-RNAz-hits.pl \
	$(top_srcdir)/configure.ac
	-$(POD2MAN) $< >$@

# locarna_rnafold_pp.1
locarna_rnafold_pp.1: $(top_srcdir)/src/locarna_rnafold_pp.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarna_rnafold_pp
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna_rnafold_pp

# locarnate.1
locarnate.1: $(top_srcdir)/Utils/locarnate $(top_srcdir)/configure.ac
	chmod +x $(top_srcdir)/Utils/locarnate
	-$(HELP2MAN) $(HELP2MANOPTS) --no-discard-stderr --output=$@ $(top_srcdir)/Utils/locarnate

# mlocarna.1
mlocarna.1: $(top_srcdir)/Utils/mlocarna $(top_srcdir)/configure.ac
	-$(POD2MAN)  $(top_srcdir)/Utils/mlocarna >$@

# mlocarna_nnames.1
mlocarna_nnames.1: $(top_srcdir)/Utils/mlocarna_nnames $(top_srcdir)/configure.ac
	-$(POD2MAN)  $(top_srcdir)/Utils/mlocarna_nnames >$@

# plot-bmprobs.1
# pp2dot.1
# reliability-profile.pl.1
reliability-profile.pl.1: $(top_srcdir)/Utils/reliability-profile.pl \
	$(top_srcdir)/configure.ac
	-$(POD2MAN) $(top_srcdir)/Utils/reliability-profile.pl >$@

# ribosum2cc.1

# sparse.1
sparse.1: $(top_srcdir)/src/sparse.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/sparse
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/sparse
