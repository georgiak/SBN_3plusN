d1:	sbnfit.cxx model.c correlation.c detector.c prob.c wrkInstance.c
	mkdir -p bkg_data
	g++ -g -std=c++11 -c prob.c -o prob.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c detector.c -o detector.o   -I $$ROOTSYS/include
	g++ -g -std=c++11 -c  sbnfit.cxx  -Wno-deprecated-declarations -I $$ROOTSYS/include
	g++ -g -std=c++11 -c  wrkInstance.c  -Wno-deprecated-declarations -I $$ROOTSYS/include
	g++ -g -std=c++11 -c model.c -o model.o -I $$ROOTSYS/include
	g++ -g -std=c++11 -c correlation.c -o correlation.o  -I $$ROOTSYS/include
	g++ -g -std=c++11  -Wall -fPIC  -Wno-deprecated-declaratios prob.o model.o detector.o correlation.o wrkInstance.o sbnfit.o -o sbnfit2 -lgsl -lgslcblas -L/scratch/markross/root/root/libMathCore.so  -lMathCore -L/scratch/markross/root/root/libMathMore.so  -lMathMore -L/scratch/markross/root/root/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics  -lThread -lMultiProc -pthread -lm -ldl -rdynamic 
