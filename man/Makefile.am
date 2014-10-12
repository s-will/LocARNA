HELP2MANOPTS=-N


dist_man_MANS=locarna.1 locarna_p.1 mlocarna.1 locarnate.1		\
         reliability-profile.pl.1 locarnap-predict-and-plot.pl.1	\
         locarnap-realign-all.pl.1 locarnap-revisit-RNAz-hits.pl.1	\
         locarnap-revcomp.pl.1 locarnap_fit.1 locarna_deviation.1	\
         exparna_p.1


locarna.1: $(top_srcdir)/src/locarna.cc
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ ../src/locarna.bin

locarna_p.1: $(top_srcdir)/src/locarna_p.cc
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ ../src/locarna_p


mlocarna.1: $(top_srcdir)/Utils/mlocarna
	-$(POD2MAN)  $(top_srcdir)/Utils/mlocarna >$@

reliability-profile.pl.1: $(top_srcdir)/Utils/reliability-profile.pl
	-$(POD2MAN) $(top_srcdir)/Utils/reliability-profile.pl >$@

locarnap-revisit-RNAz-hits.pl.1: $(top_srcdir)/Utils/locarnap-revisit-RNAz-hits.pl
	-$(POD2MAN) $< >$@

locarnap-realign-all.pl.1: $(top_srcdir)/Utils/locarnap-realign-all.pl
	-$(POD2MAN) $< >$@

locarnap-predict-and-plot.pl.1: $(top_srcdir)/Utils/locarnap-predict-and-plot.pl
	-$(POD2MAN) $< >$@

locarnap-revcomp.pl.1: $(top_srcdir)/Utils/locarnap-revcomp.pl
	-$(POD2MAN) $< >$@

locarnate.1: $(top_srcdir)/Utils/locarnate
	chmod +x $(top_srcdir)/Utils/locarnate
	-$(HELP2MAN) $(HELP2MANOPTS) --no-discard-stderr --output=$@ $(top_srcdir)/Utils/locarnate

locarna_deviation.1: $(top_srcdir)/src/deviation.cc
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ ../src/locarna_deviation

locarnap_fit.1: $(top_srcdir)/src/fit.cc
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ ../src/locarnap_fit

exparna_p.1: $(top_srcdir)/src/exparna_p.cc
	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ ../src/exparna_p

