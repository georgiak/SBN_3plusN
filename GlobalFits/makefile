d1:	src/globalFit.cxx src/ntupleProcess.cxx src/ntupleGridder.cxx src/plotter.cxx
	g++ -g -std=c++11 -c src/globalFit.cxx -o globalFit.o $(shell root-config --cflags) -I./include
	g++ -g -std=c++11 -c src/ntupleProcess.cxx -o ntupleProcess.o $(shell root-config --cflags) -I./include
	g++ -g -std=c++11 -c src/ntupleGridder.cxx -o ntupleGridder.o $(shell root-config --cflags) -I./include
	g++ -g -std=c++11 -c src/plotter.cxx -o plotter.o $(shell root-config --cflags) -I./include
	g++ -g -std=c++11  -Wall -o globalFit globalFit.o `root-config --cflags --glibs` -lMathMore -lMinuit
	g++ -g -std=c++11  -Wall -o ntupleProcess ntupleProcess.o `root-config --cflags --glibs`
	g++ -g -std=c++11  -Wall -o ntupleGridder ntupleGridder.o `root-config --cflags --glibs`
	g++ -g -std=c++11  -Wall -o plotter plotter.o `root-config --cflags --glibs`
	rm *.o
