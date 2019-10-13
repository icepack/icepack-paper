export legendre ?= .

FIGURES += $(legendre)/pressure.png

$(legendre)/pressure.png: $(legendre)/legendre.ipynb
	jupyter nbconvert $< --execute --stdout > /dev/null
