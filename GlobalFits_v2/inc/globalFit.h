#ifndef GLOBALFIT_H
#define GLOBALFIT_H

#include <sys/stat.h>
#include <unistd.h>
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph2D.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TROOT.h" //for gROOT
#include "TStyle.h" //for gStyle
#include "TH2D.h"
#include "TLegend.h"
#include "THStack.h"
#include "TImage.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TMinuit.h"
//#include "Math/Integrator.h"
//#include "Math/IntegratorMultiDim.h"
//#include "Math/Functor.h"
//#include "Math/GSLIntegrator.h"
//#include "Math/Interpolator.h"
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

/*
double sinFunc(double x){
    // Sine function for integration
    return sin(x)/x;
}
double sineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&sinFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return ig.Integral(0.,x);
}
double cosFunc(double x){
    // cosine function for integration
    return (cos(x) - 1)/x;
}
double cosineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&cosFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return TMath::EulerGamma() + log(x) + ig.Integral(0.,x);
}
*/

// probably gettin' rid of someof these fucks
struct chisqStruct{
    float chi2, chi2_det;
    void zero(){
        chi2 = 0;
        chi2_det = 0;
    }
};

struct nuVecs{
    double *EnuQE; double *NueBgr; double *Numu; double *FOsc_EnuQE, *FOsc_EnuTrue, *FOsc_LnuTrue, *FOsc_weight;
    int *NueData; int *NumuData;

    nuVecs(const int nBins,const int nBins_mu,const int nFOscEvts){
        EnuQE = new double[nBins+1];
        NueData = new int[nBins+nBins];
        NumuData = new int[nBins_mu];
        NueBgr = new double[nBins];
        Numu = new double[nBins_mu];
        FOsc_EnuQE = new double[nFOscEvts];     // reconstructed neutrino energy
        FOsc_EnuTrue = new double[nFOscEvts];   // true energy of neutrino
        FOsc_LnuTrue = new double[nFOscEvts];   // distance from production and detection points
        FOsc_weight = new double[nFOscEvts];
    }
    nuVecs(const int nBins, const int nFOscEvts){
        EnuQE = new double[nBins+1];
        NueData = new int[nBins];
        NumuData = new int[nBins];
        NueBgr = new double[nBins];
        FOsc_EnuQE = new double[nFOscEvts];     // reconstructed neutrino energy
        FOsc_EnuTrue = new double[nFOscEvts];   // true energy of neutrino
        FOsc_LnuTrue = new double[nFOscEvts];   // distance from production and detection points
        FOsc_weight = new double[nFOscEvts];
    }
};

struct boonePackage{
    std::vector < std::vector <float> > full_fractCovMatrix;
    float *EnuQE; float *NueBgr; float *Numu; float *FOsc_EnuQE, *FOsc_EnuTrue, *FOsc_LnuTrue, *FOsc_weight;
    int *NueData; int *NumuData;
	int nFOscEvts;
};
struct boonePlusPackage{
	std::vector < std::vector <float> > full_fractCovMatrix, lib_sin, lib_sinsq;
    float *NueBgr; float *Numu;
    int *NueData; int *NumuData;
	int nFOscEvts;
};
struct atmPackage{
    double * dchi2Vec, * dmuVec;
};
struct numiPackage{
    int * NueData;
    double * EnuQE, * NueBgr, * NueBgr_error;
    double * FOsc_EnuQE, * FOsc_EnuTrue, * FOsc_LnuTrue, * FOsc_weight, * FOsc_fracError;
};
struct sinSqPackage{
    std::vector < std::vector <double> > sinSqDeltaGrid;
    std::vector < std::vector <double> > sinSqDeltaGrid2;
	std::vector <double> dm2Vec;
    double * observed, * bkg, norm;
};
struct galPackage{
    double * obsRatioGal, * errorGal, * arLinesE, * arLinesBr, * arLinesXSec, * crLinesE, * crLinesBr, * crLinesXSec;
	std::vector < std::vector < double > > length, volInt;
};
struct minosncPackage{
    double theta24_bestfit;
    double sd_theta24;
};
struct minosPackage{
    double * EnuQE, * NumubarData, * NumubarBkg, * fracError, * dataErr, * EnuQE_ws, * NumubarData_ws, * NumubarBkg_ws, * fracError_ws, * dataErr_ws;
};
struct booneDisPackage{
    std::vector < std::vector <float> > full_fractCovMatrix;
    float * EnuQE, * NumuData; float *FOsc_EnuQE, *FOsc_EnuTrue, *FOsc_LnuTrue, *FOsc_weight, *LOverE;
	int nFOscEvts;
};
struct booneDisPlusPackage{
	std::vector < std::vector <float> > full_fractCovMatrix, libdis_sinsq;
	std::vector <float> libdis_noosc;
    float * NumuData;
};
struct nomadPackage{
	std::vector < std::vector <double> > sigmaRemu;
	double *observed, *bkg, *norm;
	std::vector < std::vector <double> > sinSqDeltaGrid;
    std::vector < std::vector <double> > sinSqDeltaGrid2;
};
struct ccfrPackage{
	double * observed, m_front, m_back;
	std::vector < std::vector <double> > sigmaRatio;
    std::vector < std::vector <double> > sinSqDeltaGrid_front;
    std::vector < std::vector <double> > sinSqDeltaGrid_back;
    std::vector < std::vector <double> > noOscGrid;
};
struct cdhsPackage{
	double * observed, *m_front, *m_back;
  	std::vector < std::vector <double> > sigmaRatio;
  	std::vector < std::vector <double> > sinSqDeltaGrid_front;
  	std::vector < std::vector <double> > sinSqDeltaGrid_back;
  	std::vector < std::vector <double> > noOscGrid;
	std::vector <double> dm2Vec;
};
struct bugeyPackage{
  std::vector < std::vector <double > > observed;
  std::vector < std::vector <double > > sigmaRatio;
  std::vector < std::vector <double > > energy;
  std::vector < std::vector < std::vector < double> > > sinSqDeltaGrid;

  double sigmaBigA, sigmaB, sigmaSmallA;
};
struct choozPackage{
	std::vector <double> energy;
	std::vector <double> observed;
	std::vector <double> noOsc;
	std::vector <std::vector <double > > sigmaX;

	double deltae, sigmaAlpha, sigmaG;
};
struct xsecPackage{
	std::vector <double> karmen_Enu;
	std::vector <double> karmen;
	std::vector <double> karmen_error;
	std::vector <double> lsnd_Enu;
	std::vector <double> lsnd;
	std::vector <double> lsnd_error;

	double ESigma, lsnd_sys, karmen_sys, correl_sys;
};
/*
struct minPack{
  // Here's a shitty little struct to store vars needed for minuit fit functions so I don't need to pass them in some other irritating way.
  bugeyPackage bPack;
  choozPackage cPack;
  xsecPackage xPack;
  neutrinoModel model;
};

struct integralFuncsLSND{
    double _dm2;

    double normFunc(const double * x){
        double Enu = x[1] + 1.8;
        double sigma = 0.509*sqrt(Enu);
        double michelEndpoint = 52.8;
        double xFrac = Enu/michelEndpoint;
        double xsec, flux, res;

        // Get the resolution through looking at a gaussian
        res = 1/(sqrt(2*TMath::Pi()) * sigma) * exp(-pow(x[0]-x[1],2)/(2*pow(sigma,2)));
        // Get the neutrino flux
		if(xFrac < 1) flux = pow(xFrac,2) * (3-2*xFrac);
        else flux = 0;
        // Get the neutrino cross section
        //xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.11675*pow(10,-4))*pow(Enu,2);
		xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.64283*pow(10,-4))*pow(Enu,2);
		return res * flux * xsec / pow(x[2],2);
    }
    double sinSqFunction(const double * x){
        double Enu = x[1] + 1.8;
        double sigma = 0.509*sqrt(Enu);
        double michelEndpoint = 52.8;
        double xFrac = Enu/michelEndpoint;
        double xsec, flux, res;

        // Get the resolution through looking at a gaussian
        res = 1/(sqrt(2*TMath::Pi()) * sigma) * exp(-pow(x[0]-x[1],2)/(2*pow(sigma,2)));
        // Get the neutrino flux
        if(xFrac < 1) flux = pow(xFrac,2) * (3-2*xFrac);
        else flux = 0;
        // Get the neutrino cross section
		//xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.11675*pow(10,-4))*pow(Enu,2);
		xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.64283*pow(10,-4))*pow(Enu,2);

        // Now we need to... yeah, we need to check whether oscillations are fast relative to our energy resolution and use an approximation if necessary
        double baseLineLSND = x[2];
        double EnuTrue = x[1] + 1.8;
        double delta = 1.27 * _dm2;

        double nMax = delta * baseLineLSND / (EnuTrue * TMath::Pi());
        double EnuNext = delta * baseLineLSND / ((nMax + .5) * TMath::Pi());

        if(abs(EnuTrue-EnuNext) >= 0.509 * sqrt(EnuTrue)){
			return res * flux * xsec * (1 / pow(x[2],2)) * pow(sin(delta * baseLineLSND / EnuTrue),2);
        }
        else
			return res * flux * xsec * (1 / pow(x[2],2)) * 0.5;
    }
    double sinSqFunctionCPV(const double * x){
        double Enu = x[1] + 1.8;
        double sigma = 0.509*sqrt(Enu);
        double michelEndpoint = 52.8;
        double xFrac = Enu/michelEndpoint;
        double xsec, flux, res;

        // Get the resolution through looking at a gaussian
        res = 1/(sqrt(2*TMath::Pi()) * sigma) * exp(-pow(x[0]-x[1],2)/(2*pow(sigma,2)));
        // Get the neutrino flux
        if(xFrac < 1) flux = pow(xFrac,2) * (3-2*xFrac);
        else flux = 0;
        // Get the neutrino cross section
		//xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.11675*pow(10,-4))*pow(Enu,2);
		xsec = (-0.42222*pow(10,-3)) + (0.64283*pow(10,-4))*Enu + (0.64283*pow(10,-4))*pow(Enu,2);

        // Now we need to... yeah, we need to check whether oscillations are fast relative to our energy resolution and use an approximation if necessary
        double baseLineLSND = x[2];
        double EnuTrue = x[1] + 1.8;
        double delta = 1.27 * _dm2;

        double nMax = delta * baseLineLSND / (EnuTrue * TMath::Pi());
        double EnuNext = delta * baseLineLSND / ((nMax + .5) * TMath::Pi());
		if(abs(EnuTrue-EnuNext) >= 0.509 * sqrt(EnuTrue)){
            return res * flux * xsec * (1 / pow(x[2],2)) * sin(2 * delta * baseLineLSND / EnuTrue);
        }
        else    return 0;
    }
};

struct integralFuncsNomad{
    double _dm2, _EnuAvg;

    double sinSqFunction(const double x){
        double l1 = .422; double l2 = .837; double baseline = (l1 + l2)/2;

        double a = 1.27 * _dm2;
        double arg1 = a * l1 / x;
        double arg2 = a * l2 / x;
        double num = (TMath::Pi() / 2.) * (1. / (1.27 * _dm2 * baseline)) * pow(_EnuAvg,2);
        double den = 1. + (TMath::Pi() / 2.) * (1. / (1.27 * _dm2 * baseline)) * _EnuAvg;
        double deltaEnu = num / den;
        // Check if oscillations are very fast relative to energy resolution
        if(deltaEnu >= .032 * sqrt(_EnuAvg) + .01 * _EnuAvg){
            double term;
            if(a > 0)
                term = .5 * (1. - .5 * (sin(2. * arg1) - sin(2. * arg2)) / (arg1 - arg2)) * (l2 - l1);
            else term = 0.;
            return term;
        }
       else return .5 * (l2 - l1);
   }
    double sinSqFunctionCPV(const double x){
        double l1 = .422; double l2 = .837; double baseline = (l1 + l2)/2;

        double a = 1.27 * _dm2;
        double arg1 = a * l1 / x;
        double arg2 = a * l2 / x;
        double num = (TMath::Pi() / 2.) * (1. / (1.27 * _dm2 * baseline)) * pow(_EnuAvg,2);
        double den = 1. + (TMath::Pi() / 2.) * (1. / (1.27 * _dm2 * baseline)) * _EnuAvg;
        double deltaEnu = num / den;
        // Check if oscillations are very fast relative to energy resolution
        if(deltaEnu >= .032 * sqrt(_EnuAvg) + .01 * _EnuAvg){
            double term;
            if(a > 0)
                term = x / (2. * a) * (cos(2. * arg1) - cos(2. * arg2));
            else term = 0.;
            return term;
        }
        else return 0.;
    }
};
struct integralFuncsCDHS1{
    // This could have been dealt with better, but it wasn't.
    bool _noOsc, _front;
    double _dm2, _emu;

    double intFunc1(const double x){

        double lBack = 0.885; double lFront = 0.130;
        double deltaBack = .144; double deltaFront = .074;
        double l1, l2;
        if(!_front){
            l1 = lBack - deltaBack / 2.;
            l2 = lBack + deltaBack / 2.;
        }
        else{
            l1 = lFront - deltaFront / 2.;
            l2 = lFront + deltaFront / 2.;
        }

        double noOsc_term = (l2 - l1) / (l1 * l2);
        double delta = 1.27 * _dm2;// / x;
        double flux = pow(10,6) * exp(-x/1.);
        double xsec;
        double fitPar[3] = {0.97443, 0.72555, 0.024532};
        if(_emu <= fitPar[0] * x)   xsec = TMath::Max(exp(-pow(_emu - fitPar[0]*x,2)/(2.*pow(x*fitPar[1],2))),0.);
        else    xsec = TMath::Max(exp(-pow(_emu - fitPar[0]*x,2)/(2.*pow(x*fitPar[2],2))),0.);

        if(!_noOsc){
			double num1 = x - x * cos(2. * delta * l1 / x) - 2. * delta * l1 * sineInt(2. * delta * l1 / x);
            double num2 = x - x * cos(2. * delta * l2 / x) - 2. * delta * l2 * sineInt(2. * delta * l2 / x);
            double den1 = 2. * x * l1;  double den2 = 2. * x * l2;
            double term1 = num1 / den1;    double term2 = num2 / den2;

			//std::cout << "TERM1-TERM2: " << term1-term2 << std::endl;
            return (term1 - term2) * flux * xsec;
        }
        else return noOsc_term * flux * xsec;
    }
};
struct integralFuncsCDHS2{
    bool _noOsc, _front;
    double _dm2;

    double intFunc2(const double x){
        integralFuncsCDHS1 integralFuncs;
        integralFuncs._noOsc = _noOsc;
        integralFuncs._emu = x;
        if(!_noOsc) integralFuncs._dm2 = _dm2;

        ROOT::Math::Functor1D wf(&integralFuncs,&integralFuncsCDHS1::intFunc1);
        ROOT::Math::Integrator ig;
        ig.SetFunction(wf);

        return ig.Integral(x, 7.);
    }
};
struct integralFuncsBugey{
  int _jB;
  double _dm2, _energy;

  double sinSqFunction(const double x){

      double sigma = .4;
      double baselines[] = {15., 40., 95.};
      double EMean = _energy;
      double l1 = baselines[_jB] - 1.5;
      double l2 = baselines[_jB] + 1.5;
      double a = 1.27 * _dm2;
      double resolFunction = (1./sqrt(2*TMath::Pi()*pow(sigma,2))) * exp(-pow(x - EMean,2)/(2*pow(sigma,2)));
      double num1 = x - x * cos(2. * a * l1 / x) - 2. * a * l1 * sineInt(2. * a * l1 / x);
      double num2 = x - x * cos(2. * a * l2 / x) - 2. * a * l2 * sineInt(2. * a * l2 / x);
      double den1 = 2. * x * l1;  double den2 = 2. * x * l2;
      double term1 = num1 / den1; double term2 = num2 / den2;

      return resolFunction * (term1 - term2) * l1 * l2 / (l2 - l1);
  }
};
*/

bool jobOpt();

#endif
