export legendre ?= .

FIGURES += $(legendre)/pressure.png

$(legendre)/pressure.png: $(legendre)/legendre.ipynb
	jupyter nbconvert $< --to ipynb --execute --stdout > /dev/null
