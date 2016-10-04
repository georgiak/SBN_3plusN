<<<<<<< HEAD
d1:	sbnfit.cxx model.c correlation.c detector.c prob.h
	mkdir -p bkg_data
	g++ -g -std=c++11 -c prob.c -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c detector.c -o detector.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c sbnfit.cxx -o sbnfit.o  -I $$ROOTSYS/include -Wno-deprecated-declarations
	g++ -g -std=c++11 -c model.c -o model.o -I $$ROOTSYS/include
	g++ -g -std=c++11 -c correlation.c -o correlation.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -Wall -fPIC  -Wno-deprecated-declarations -I $$ROOTSYS/include -o sbnfit sbnfit.o prob.o model.o detector.o correlation.o `root-config --libs`  -lgsl -lgslcblas 
=======
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
>>>>>>> 761985f9c500b2ffc6cf62deae039608478c9459
