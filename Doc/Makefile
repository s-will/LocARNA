FIGS=Figs/relprof.pdf

all:	locarna-doc.pdf

%.output: %.fa
	mlocarna --probabilistic --consistency-transformation $< >$@

example.evaluation: example.output
	mlocarna --evaluate example.out/results/result.aln example.fa  >$@

example.fa:
	cp ../Examples/archaea.fa example.fa

example2.fa:
	cp ../Examples/tRNA_5.fa example2.fa

locarna-doc.pdf: locarna-doc.tex $(FIGS) example.fa example.output
	pdflatex $<
	pdflatex $<

Figs/relprof.pdf: example.output
	reliability-profile.pl example.out --out=Figs/relprof.pdf
