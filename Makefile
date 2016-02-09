d1:	sbnfit.cxx model.c correlation.c 
	g++ -g -std=c++11 -c sbnfit.cxx -o sbnfit.o -I/usr/include/root 
	g++ -g -std=c++11 -c model.c -o model.o -I. -I/usr/include/root  
	g++ -g -std=c++11 -c correlation.c -o correlation.o -I. -I/usr/include/root  
	g++ -g -std=c++11 -o sbnfit sbnfit.o model.o correlation.o -lgomp -lnlopt -lgsl -lgslcblas  -L/usr/lib/x86_64-linux-gnu -lGui -lCore -lCint -lRIO -lNet -lHist  -lTree -lRint  -lMatrix -lMathCore -lThread -pthread -lm -ldl -rdynamic 

	rm *.o

