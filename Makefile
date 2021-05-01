
.PHONY: clean

all: icepack.pdf

src ?= .
legendre := $(src)/demos/legendre
sliding := $(src)/demos/sliding
mismip := $(src)/demos/mismip
larsen := $(src)/demos/larsen
gibbous := $(src)/demos/gibbous
icesheet := $(src)/demos/ice-sheet

include $(legendre)/Makefile
include $(sliding)/Makefile
include $(mismip)/Makefile
include $(larsen)/Makefile
include $(gibbous)/Makefile
include $(icesheet)/Makefile

icepack.pdf: icepack.tex icepack.bib $(FIGURES)
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

INITIAL_SUBMISSION_COMMIT=3554c618468320c06c25bca46fd2f97c5d1e860c

icepack-initial-submission.bib:
	git show $(INITIAL_SUBMISSION_COMMIT):icepack.bib > $@

icepack-initial-submission.tex: icepack-initial-submission.bib
	git show $(INITIAL_SUBMISSION_COMMIT):icepack.tex > $@
	sed -i s/icepack\.bib/icepack-initial-submission\.bib/g $@

icepack-initial-submission.pdf: icepack-initial-submission.tex icepack-initial-submission.bib $(FIGURES)
	pdflatex icepack-initial-submission
	bibtex icepack-initial-submission
	pdflatex icepack-initial-submission
	pdflatex icepack-initial-submission

LATEXDIFF_OPTS=--type CCHANGEBAR --math-markup=off

icepack-diff.tex: icepack.tex icepack-initial-submission.tex
	latexdiff $(LATEXDIFF_OPTS) icepack-initial-submission.tex icepack.tex > icepack-diff.tex
	latexdiff $(LATEXDIFF_OPTS) icepack-initial-submission.bbl icepack.bbl > icepack-diff.bbl

icepack-diff.pdf: icepack-diff.tex
	pdflatex icepack-diff
	bibtex icepack-diff
	pdflatex icepack-diff
	pdflatex icepack-diff

clean:
	rm *.pdf $(FIGURES)
