
.PHONY: clean

all: icepack.pdf

demos/legendre/pressure.png: demos/legendre/legendre.ipynb
	jupyter nbconvert $< --execute --stdout > /dev/null

demos/mismip/steady-state-level0.h5: demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --time 6e3 --timestep 1.0 --verbose

demos/mismip/steady-state-level1.h5: demos/mismip/steady-state-level0.h5 demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --output-level 1 --input $< --input-level 0 --time 3.6e3 --timestep 1.0 --verbose

demos/mismip/steady-state-level2.h5: demos/mismip/steady-state-level1.h5 demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --output-level 2 --input $< --input-level 1 --time 1000 --timestep 1.0 --verbose

demos/mismip/steady-state.h5: demos/mismip/steady-state-level2.h5 demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --output-level 3 --input $< --input-level 2 --time 200 --timestep 1.0 --verbose

demos/mismip/retreated.h5: demos/mismip/steady-state.h5 demos/mismip/mismip+.py
	python3 demos/mismip/mismip+.py --output $@ --output-level 3 --input $< --input-level 3 --time 100 --timestep 0.04166666666666666 --melt --verbose

demos/sliding/sliding-law.png: demos/sliding/sliding-law.py
	python3 demos/sliding/sliding-law.py --output $@

demos/mismip/%.png: demos/mismip/%.h5 demos/mismip/plot.py
	python3 demos/mismip/plot.py --input $< --level 3 --output $@

FIGURES=demos/legendre/pressure.png demos/sliding/sliding-law.png demos/mismip/steady-state.png demos/mismip/retreated.png

icepack.pdf: icepack.tex icepack.bib $(FIGURES)
	pdflatex icepack
	bibtex icepack
	pdflatex icepack
	pdflatex icepack

clean:
	rm *.pdf $(FIGURES)
