d1:	sbnfit.cxx model.c correlation.c detector.c prob.h
	mkdir -p bkg_data
	g++ -g -std=c++11 -c prob.c -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c detector.c -o detector.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c sbnfit.cxx -o sbnfit.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -c model.c -o model.o -I $$ROOTSYS/include
	g++ -g -std=c++11 -c correlation.c -o correlation.o  -I $$ROOTSYS/include
	g++ -g -std=c++11 -Wall -fPIC   -I $$ROOTSYS/include -o sbnfit sbnfit.o prob.o model.o detector.o correlation.o `root-config --libs`  -lgsl -lgslcblas 
