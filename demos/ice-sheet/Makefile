export icesheet ?= .

SCRIPT := $(icesheet)/ice-sheet.py

FIGURES += $(icesheet)/ice-sheet.png

$(icesheet)/ice-sheet.h5: $(SCRIPT)
	python3 $(SCRIPT) --timestep 10 --num-steps 50 --output $@

$(icesheet)/ice-sheet.png: $(icesheet)/ice-sheet.h5 $(icesheet)/make_plots.py
	python3 $(icesheet)/make_plots.py --input $< --output $@
