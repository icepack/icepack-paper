
.PHONY: clean

all: icepack.pdf

demos/legendre/pressure.png: demos/legendre/legendre.ipynb
	jupyter nbconvert $< --execute --stdout > /dev/null

FIGURES=demos/legendre/pressure.png

icepack.pdf: icepack.tex icepack.bib $(FIGURES)
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

clean:
	rm *.pdf $(FIGURES)
