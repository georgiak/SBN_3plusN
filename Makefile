SHELL := /bin/bash

MB:	main.tex
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex
	rm *.log
	rm main.aux
	rm main.out
	rm main.toc
	rm main.bbl
	rm main.blg
	mv main.pdf SBN.pdf	
