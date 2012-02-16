RFILES=$(wildcard *.R)
SKFILES=$(wildcard *.sk)
SVGFILES=$(wildcard *.svg)
SVGZFILES=$(wildcard *.svgz)
PNGFILES=$(wildcard *.png)
JPGFILES=$(wildcard *.jpg)

DPI=150 ## dpi for conversion svg -> pdf

OVERLAYPDFS=

TGTPDFFILES=$(SKFILES:.sk=.pdf) $(SVGFILES:.svg=.pdf) \
	 $(PNGFILES:.png=.pdf) \
	 $(OVERLAYPDFS) \
	 $(JPGFILES:.jpg=.pdf) \
	 $(RFILES:.R=.pdf) $(SVGZFILES:.svgz=.pdf) \

all: $(TGTPDFFILES)

$.pdf: %.jpg
	convert $< $@

%.pdf: %.R
	R --vanilla <$<

%.eps: %.sk
	sk2ps $< >$@

%.pdf: %.eps
	epstopdf $<

%.pdf: %.png
	convert $< $@

%.pdf: %.jpg
	convert $< $@

%.pdf: %.png
	convert $< $@

%.pdf: %.jpg
	convert $< $@

%.pdf: %.svgz
	inkscape $< --export-pdf=$@

%.png: %.svgz
	inkscape $< --export-dpi=$(DPI) --export-png=$@

%.pdf: %.svg
	inkscape $< --export-pdf=$@

%.png: %.svg
	inkscape $< --export-dpi=$(DPI) --export-png=$@


clean:
	$(RM) $(TGTPDFFILES) $(OVERLAYPDFS:.pdf=.svg)
