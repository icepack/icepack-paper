export mismip ?= .

FIGURES += $(mismip)/steady-state.png $(mismip)/retreated.png

SCRIPT := $(mismip)/mismip+.py

$(mismip)/steady-state-level0.h5: $(SCRIPT)
	python3 $(SCRIPT) --output $@ --time 6e3 --timestep 1.0 --verbose

$(mismip)/steady-state-level1.h5: $(mismip)/steady-state-level0.h5 $(SCRIPT)
	python3 $(SCRIPT) --verbose \
	    --output $@ --output-level 1 \
	    --input $< --input-level 0 \
	    --time 3.6e3 --timestep 1.0

$(mismip)/steady-state-level2.h5: $(mismip)/steady-state-level1.h5 $(SCRIPT)
	python3 $(SCRIPT) --verbose \
	    --output $@ --output-level 2 \
	    --input $< --input-level 1 \
	    --time 1e3 --timestep 1.0

$(mismip)/steady-state.h5: $(mismip)/steady-state-level2.h5 $(SCRIPT)
	python3 $(SCRIPT) --verbose \
	    --output $@ --output-level 3 \
	    --input $< --input-level 2 \
	    --time 200 --timestep 1.0

$(mismip)/retreated.h5: $(mismip)/steady-state.h5 $(SCRIPT)
	python3 $(SCRIPT) --verbose \
	    --output $@ --output-level 3 \
	    --input $< --input-level 3 \
	    --time 100 0.04166666666666666 --melt

$(mismip)/%.png: $(mismip)/%.h5 $(mismip)/plot.py
	python3 $(mismip)/plot.py --input $< --level 3 --output $@
