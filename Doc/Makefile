FIGS=Figs/relprof.pdf

MANS=$(wildcard ../man/*)
DOCPDFS=$(subst ../,,$(MANS:.1=.pdf))

.PHONY: manpagepdfs

all:	locarnap-doc.pdf

%.output: %.fa
	mlocarna --probabilistic --consistency-transformation $< >$@

example.evaluation: example.output
	mlocarna --evaluate example.output/results/result.aln example.fa  >$@

example.fa:
	cp ../Examples/archaea.fa example.fa

example2.fa:
	cp ../Examples/tRNA_5.fa example2.fa

locarnap-doc.pdf: locarnap-doc.tex $(FIGS) example.fa example.output manpagepdfs
	pdflatex $<
	pdflatex $<
# append all manual pages
	pdfjoin $@ $(DOCPDFS) --outfile $@-joined
	mv $@-joined $@

Figs/relprof.pdf: example.output
	reliability-profile.pl example.output --out=Figs/relprof.pdf


manpagepdfs: $(DOCPDFS)
	if [ ! -d man ] ; then mkdir man ; fi

man/%.pdf: ../man/%.1
	man -t $< | ps2pdf - > $@


clean:
	$(RM) locarnap-doc.aux locarnap-doc.log locarnap-doc.out locarnap-doc.pdf
