
.PHONY: clean

all: icepack.pdf

demos/legendre/pressure.png: demos/legendre/legendre.ipynb
	jupyter nbconvert $< --execute --stdout > /dev/null

demos/mismip/steady-state.h5: demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --time 10e3 --timestep 1.0 --verbose | tee log.txt

demos/mismip/retreated.h5: demos/mismip/steady-state.h5 demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --input $< --time 100 --timestep 0.04166666666666666 --melt --verbose | tee log.txt

demos/sliding/sliding-law.png: demos/sliding/sliding-law.py
	python3 demos/sliding/sliding-law.py --output $@

demos/mismip/%.png: demos/mismip/%.h5 demos/mismip/plot.py
	python3 demos/mismip/plot.py --input $< --output $@

FIGURES=demos/legendre/pressure.png demos/sliding/sliding-law.png demos/mismip/steady-state.png demos/mismip/retreated.png

icepack.pdf: icepack.tex icepack.bib $(FIGURES)
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

clean:
	rm *.pdf $(FIGURES)
