HELP2MANOPTS=-Nl

## man pages for binaries are generated using help2man 
## man pages for perl scripts are generated using pod2man 

## generation of man pages is adapted from the gnu hello example:
## generate man pages only in vc-controlled dirs, i.e. avoid man page
## generation when compiling from distributed source (no pages are
## generated when cross compiling; this is overly strict for perl
## scripts, but we don't really care for now)


EXTRA_DIST =
CLEANFILES =

## Binaries, where man pages shall be generated using help2man
##
## auto-generate list by
## file * | grep 'ELF' | cut -d: -f1 | tr '\n' ' '
## in (locarna only) installation bin directory
##
## take out locarna for special handling (due to naming fix)
##
help2man_prgs1=exparna_p locarna_deviation locarna_p locarnap_fit	\
	locarna_rnafold_pp ribosum2cc sparse
help2man_prgs=$(help2man_prgs1) locarna

## Perl scripts, where man pages shall be generated using pod2man
##
## auto-generate list by
## file * | grep 'Perl script' | cut -d: -f1 | tr '\n' ' '
## in (locarna only) installation bin directory
##
pod2man_prgs=aln2fa.pl alnsel.pl aln-seqs.pl average-dot.pl dot2pp	\
exploc_p gen-reliab-dot.pl locarna-mea.pl locarna-motif-scan		\
locarnap-predict-and-plot.pl locarnap-realign-all.pl			\
locarnap-revcomp.pl locarnap-revisit-RNAz-hits.pl locarnate mlocarna	\
mlocarna_nnames plot-bmprobs pp2dot reliability-profile.pl

if !VC_CROSS_COMPILING

man_MANS=$(help2man_prgs:=.1) $(pod2man_prgs:=.1)

EXTRA_DIST += $(man_MANS)
endif

if BUILD_FROM_VC

## For portability, we generate single rules for each man page.  Using
## non-portable gnu make-specific pattern rules, the same effect would
## be achieved by
# %.1:	$(top_builddir)/src/%
# 	-$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $<

# %.1: $(top_srcdir)/Utils/%
# 	-$(POD2MAN) $< >$@

include help2man-rules pod2man-rules
CLEANFILES += help2man-rules pod2man-rules

help2man-rules:
	for i in $(help2man_prgs1); do\
	  printf "$$i.1:	$(top_builddir)/src/$$i\n\
	\t$(HELP2MAN) $(HELP2MANOPTS) --output=\$$@ $$<\n\n" ;\
	done > $@

pod2man-rules:
	for i in $(pod2man_prgs); do\
	  printf "$$i.1:	$(top_srcdir)/Utils/$$i\n\
	\t$(POD2MAN) $$< \$$@\n\n" ;\
	done > $@

## overwrite the rule for locarna
locarna.1: $(top_builddir)/src/locarna.bin
	$(HELP2MAN) $(HELP2MANOPTS) --output=$@ $(top_builddir)/src/locarna.bin
#	strip .bin (special name handling for locarna binary)
	sed -i 's/locarna.bin/locarna/g' $@

CLEANFILES += $(man_MANS)
endif
