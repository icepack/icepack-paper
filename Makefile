
.PHONY: clean

all: icepack.pdf

icepack.pdf: icepack.tex icepack.bib
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

clean:
	rm *.pdf
