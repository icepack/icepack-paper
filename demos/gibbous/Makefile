export gibbous ?= .

FIGURES += $(gibbous)/damage.png

SCRIPT := $(gibbous)/ice-shelf.py

$(gibbous)/ice-shelf.geo: $(gibbous)/make_mesh.py
	python3 $< --output $@

$(gibbous)/ice-shelf.msh: $(gibbous)/ice-shelf.geo
	gmsh -2 -format msh2 -v 0 -o $@ $<

$(gibbous)/steady-state-undamaged.h5: $(gibbous)/ice-shelf.msh $(SCRIPT)
	python3 $(SCRIPT) --mesh $< --output $@ --final-time 400 --num-steps 200

$(gibbous)/steady-state-damaged.h5: $(gibbous)/steady-state-undamaged.h5 $(SCRIPT)
	python3 $(SCRIPT) --mesh $(gibbous)/ice-shelf.msh --input $< --output $@ \
	    --damage --final-time 400 --num-steps 200

$(gibbous)/damage.png: $(gibbous)/steady-state-undamaged.h5 $(gibbous)/steady-state-damaged.h5 $(gibbous)/make_plots.py
	python3 $(gibbous)/make_plots.py --mesh $(gibbous)/ice-shelf.msh --undamaged $(gibbous)/steady-state-undamaged --damaged $(gibbous)/steady-state-damaged --output $@
