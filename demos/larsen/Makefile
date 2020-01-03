export larsen ?= .

FIGURES += $(larsen)/parameter.png

$(larsen)/parameter.h5: $(larsen)/larsen.py
	python3 $(larsen)/larsen.py \
	    --geometry $(larsen)/larsen.geo \
	    --mesh $(larsen)/larsen.msh --level 1 \
	    --thickness-smoothing 8e3 --velocity-smoothing 2.5e3 \
	    --regularization 7.5e2 --gamma 1.0 \
	    --output $@

$(larsen)/%.png: $(larsen)/parameter.h5 $(larsen)/plot.py
	python3 $(larsen)/plot.py \
	    --mesh $(larsen)/larsen.msh --level 1 \
	    --input $(larsen)/parameter.h5 \
	    --output $(larsen)/parameter.png
