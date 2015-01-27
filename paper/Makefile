all: Diversity_MS.pdf

Diversity_MS.pdf: Diversity_MS.tex
	pdflatex Diversity_MS
	bibtex Diversity_MS
	pdflatex Diversity_MS
	pdflatex Diversity_MS

tidy:
	rm -f Diversity_MS.{aux,bbl,blg,log}

clean: tidy
	rm -f Diversity_MS.pdf

.PHONY: all tidy clean
