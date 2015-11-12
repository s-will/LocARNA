HELP2MANOPTS=-Nl

EXTRA_DIST =
CLEANFILES =

if !VC_CROSS_COMPILING
dist_man_MANS=aln2fa.pl.1 alnsel.pl.1 aln-seqs.pl.1 average-dot.pl.1	\
	dot2pp.1 exparna_p.1 exploc_p.1 gen-reliab-dot.pl.1 locarna.1	\
	locarna_deviation.1 locarna-mea.pl.1 locarna-motif-scan.1	\
	locarna_p.1 locarnap_fit.1 locarna.pl.1				\
	locarnap-predict-and-plot.pl.1 locarnap-realign-all.pl.1	\
	locarnap-revcomp.pl.1 locarnap-revisit-RNAz-hits.pl.1		\
	locarna_rnafold_pp.1 locarnate.1 mlocarna.1 mlocarna_nnames.1	\
	plot-bmprobs.1 pp2dot.1 reliability-profile.pl.1 sparse.1	\
	ribosum2cc.1

EXTRA_DIST += $(man_MANS)
endif

if BUILD_FROM_VC

%.1:	$(top_builddir)/src/%
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $<

%.1: $(top_srcdir)/Utils/%
	-$(POD2MAN) $< >$@


locarna.1: $(top_srcdir)/src/locarna.cc $(top_srcdir)/configure.ac \
	| $(top_builddir)/src/locarna.bin
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna.bin
#	strip .bin (special name handling for locarna binary)
	sed -i 's/locarna.bin/locarna/g' $@

CLEANFILES += $(man_MANS)
endif
