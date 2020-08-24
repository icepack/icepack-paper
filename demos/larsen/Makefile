export larsen ?= .

FIGURES += $(larsen)/parameter.png

$(larsen)/parameter.h5: $(larsen)/inverse.py
	python3 $(larsen)/inverse.py \
	    --geometry $(larsen)/larsen.geo \
	    --mesh $(larsen)/larsen.msh \
	    --thickness-smoothing 2e3 \
	    --regularization 7.5e3 \
	    --output $@

$(larsen)/%.png: $(larsen)/parameter.h5 $(larsen)/plot.py
	python3 $(larsen)/plot.py \
	    --mesh $(larsen)/larsen.msh \
	    --input $(larsen)/parameter.h5 \
	    --output $(larsen)/parameter.png
