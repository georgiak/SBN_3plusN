#ifndef WRK_H_
#define WRK_H_


#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TError.h"

#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"

#include "params.h"
#include "prob.h"
#include "model.h"
#include "correlation.h"
#include "wrkInstance.h"


class wrkInstance {
	
		int matrix_size ;
		int matrix_size_c ;
		int bigMsize;
		int contMsize;

	public:



	int beam_mode; // 0 is nu only 1 is nubar+nu
	int which_mode; //app, dis or both
	bool isVerbose;
	double pot;
	double pot_bar;

	neutrinoModel nullModel;
	neutrinoModel workingModel;

/*	SBN_spectrum * bkgspec;
	SBN_spectrum * SigSpec;

	SBN_spectrum * bkgbarspec;
	SBN_spectrum * SigBarSpec;
*/
	std::vector<double > back6 ;
	std::vector<double > back  ;

	std::vector<double > backbar6 ;
	std::vector<double > backbar  ;

	std::vector<double > pred6;
	std::vector<double > predbar6;
	std::vector<double >  pred_all_12;

	std::vector<double> back_all;
	std::vector<double> back_all_12;

	double Current_Chi;

	std::vector<std::vector<double >> vMcI;
	std::vector<std::vector<double >> vMc;

	ROOT::Math::Minimizer* min ;     

/*   ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );
 
   min.SetMaxFunctionCalls(100000);
   min.SetMaxIterations(10000);
   min.SetTolerance(0.01);
 
   ROOT::Math::Functor f(&RosenBrock,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = { -1.,1.2};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   //    min.SetVariable(0,"x",variable[0], step[0]);
   //       min.SetVariable(1,"y",variable[1], step[1]);
   //        
   //           min.Minimize(); 
   //            
   //               const double *xs = min.X();
   //                  cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
   //                          << RosenBrock(xs) << endl;
   //                           
   //                              return 0;
*/



	wrkInstance(int channel_mode, int fbeam_mode , double pot_scale, double pot_scale_bar); //for pot analysis
	wrkInstance(int channel_mode, int fbeam_mode);
	~wrkInstance();

	
	wrkInstance(neutrinoModel signalModel, int channel_mode, int fbeam_mode , double pot_scale, double pot_scale_bar);
	double inject_signal(neutrinoModel signalModel, int channel_mode, int fbeam_mode , double pot_scale, double pot_scale_bar);

	double calc_chi(neutrinoModel signalModel, int runnumber);
	double calc_chi(neutrinoModel signalModel, int runnumber, double pot_scale, double pot_scale_bar);
	double minim_calc_chi(const double * xx);

	double calc_chi_POT_vector(neutrinoModel newModel, std::vector<double> vecin , int runnumber, double potin, double potinbar);

	int init_minim();
	double minimize(double phi45, double ipot, double ipotbar);
	double minimize(neutrinoModel newModel, double ipot, double ipotbar);
	int reset_minim();

	int clear_all();
};
#endif
