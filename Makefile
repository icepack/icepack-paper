
.PHONY: clean

all: icepack.pdf

src ?= .
legendre := $(src)/demos/legendre
sliding := $(src)/demos/sliding
mismip := $(src)/demos/mismip
larsen := $(src)/demos/larsen
gibbous := $(src)/demos/gibbous

include $(legendre)/Makefile
include $(sliding)/Makefile
include $(mismip)/Makefile
include $(larsen)/Makefile
include $(gibbous)/Makefile

icepack.pdf: icepack.tex icepack.bib $(FIGURES)
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

clean:
	rm *.pdf $(FIGURES)
