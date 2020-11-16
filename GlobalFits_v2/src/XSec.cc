#include "XSec.h"

neutrinoModel minModelXSec;
double ESigma, Lsnd_sys, Karmen_sys, Correl_sys;
std::vector < double > Karmen_Enu, Karmen, Karmen_error, Lsnd_Enu, Lsnd, Lsnd_error;

int XSec::Init(std::string dataLoc, Oscillator osc, bool debug){

  Karmen_Enu.push_back(28.70);
  Karmen_Enu.push_back(32.71);
  Karmen_Enu.push_back(36.55);
  Karmen_Enu.push_back(40.90);
  Karmen_Enu.push_back(45.42);
  Karmen_Enu.push_back(49.51);

  Karmen.push_back(4.68);
  Karmen.push_back(6.31);
  Karmen.push_back(9.18);
  Karmen.push_back(13.17);
  Karmen.push_back(23.57);
  Karmen.push_back(36.77);

  Karmen_error.push_back(1.26);
  Karmen_error.push_back(1.39);
  Karmen_error.push_back(1.52);
  Karmen_error.push_back(1.83);
  Karmen_error.push_back(2.47);
  Karmen_error.push_back(5.76);

  Lsnd_Enu.push_back(37.80);
  Lsnd_Enu.push_back(41.48);
  Lsnd_Enu.push_back(44.93);
  Lsnd_Enu.push_back(47.89);
  Lsnd_Enu.push_back(49.91);

  Lsnd.push_back(11.42);
  Lsnd.push_back(17.59);
  Lsnd.push_back(21.39);
  Lsnd.push_back(25.20);
  Lsnd.push_back(35.17);

  Lsnd_error.push_back(0.92);
  Lsnd_error.push_back(1.31);
  Lsnd_error.push_back(1.58);
  Lsnd_error.push_back(2.23);
  Lsnd_error.push_back(4.99);

  ESigma = 0.15;

  Lsnd_sys = sqrt(pow(.1,2) - pow(.07,2));
  Karmen_sys = sqrt(pow(.088,2) - pow(.07,2));
  Correl_sys = sqrt(pow(.12,2) + pow(.07,2));

  // Initialize Minuit minimizer
  gMinuit = new TMinuit(6);

  dof = 11;

  // Iniitialize Output tree
  chi2Nt = new OutTree("XSec");

  if(debug) std::cout << "XSec initialized. Bins: " << 11 << std::endl;

  return dof;
}


float XSec::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  double chisq = 0;
  minModelXSec = model;
  double kkarmen, klsnd, kcorrel;

  double arglis[2];
  int ierflag;
  arglis[0] = -1.;
  gMinuit->SetFCN(fcnXSec);
  // Okay, let's get this minuit garbage started
  gMinuit->mnexcm("SET PRI",arglis,1,ierflag);
  gMinuit->mnexcm("SET NOWARNING",arglis,0,ierflag);
  arglis[0] = 1.;
  gMinuit->mnexcm("SET STR",arglis,1,ierflag);
  // Now, clear parameters and set 'em anew
  gMinuit->mnexcm("CLE",arglis,1,ierflag);
  gMinuit->mnparm(0,TString("k_lsnd"),1.,.1,0.,0.,ierflag);
  gMinuit->mnparm(1,TString("k_karmen"),1.,.1,0.,0.,ierflag);
  gMinuit->mnparm(2,TString("k_correl"),1.,.1,0.,0.,ierflag);

  // With these high dm2's, CHOOZ oscillations average out so we don't need to fit energy scale factor g
  arglis[0] = 5.;
  gMinuit->mnexcm("CALL", arglis,1,ierflag);
  double disttomin, errdef;   int npari, nparx, istat;
  arglis[0] = 100000.;
  arglis[1] = .1;
  gMinuit->mnexcm("MINI", arglis,1,ierflag);

  gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
  chi2 = chisq;

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug) std::cout << "XSec Chi2: " << chi2 << std::endl;

  return chi2;
}

// Important functions for chi2 calculation
double getSinSqTerm(double dm2, double E1, double E2, double l1, double l2){

	double k = 1.27 * dm2;
	double a11 = 2. * k * l1/E1;
	double a12 = 2. * k * l1/E2;
	double a21 = 2. * k * l2/E1;
	double a22 = 2. * k * l2/E2;

	double sinSqTerm = -4 * E2 * k * l2 - 4 * pow(k,2) * sineInt(a21) * pow(l2,2)
			+ 4. * pow(k,2) * sineInt(a22) * pow(l2,2)
			+ 4. * E1 * k * l2 - sin(a21) * pow(E1,2) + sin(a22) * pow(E2,2)
			- 2. * k * cos(a21) * E1 * l2 + 2. * k * cos(a22) * E2 * l2
			+ 4. * E2 * k * l1 + 4. * pow(k,2) * sineInt(a11) * pow(l1,2)
			- 4. * pow(k,2) * sineInt(a12) * pow(l1,2) - 4. * E1 * k * l1 + sin(a11) * pow(E1,2)
			- sin(a12) * pow(E2,2) + 2. * k * cos(a11) * E1 * l1 - 2. * k * cos(a12) * E2 * l1;

	sinSqTerm = -sinSqTerm / (8. * k) / (E2 - E1) / (l2 - l1);

	// Now, check energy resolution and return .5 if necessary
	double EAvg = (E1 + E2)/2.;	double lAvg = (l1 + l2)/2.;
	double n = 1.27 * dm2 * lAvg / (EAvg * TMath::Pi());
	double ENext = 1.27 * dm2 * lAvg / ((n + .5) * TMath::Pi());
	if(n > 2 && abs(EAvg - ENext)/EAvg < ESigma)
		sinSqTerm = 0.5;

	return sinSqTerm;
}
double osc_int(double E1, double E2, double l1, double l2){

  double sin22th = minModelXSec.ProbAmp("ee");
  double dm2 = minModelXSec.Dm2();

  double osc_int = sin22th * getSinSqTerm(dm2,E1,E2,l1,l2);
  return osc_int;

	/*  Drastically simplify to 3+1  model
	if(minModelXSec.dm41Sq != 0)
		sin41sq = getSinSqTerm(minModelXSec.dm41Sq,E1,E2,l1,l2);
	else
		sin41sq = 0;

	if(minModelXSec.dm51Sq != 0){
		sin51sq = getSinSqTerm(minModelXSec.dm51Sq,E1,E2,l1,l2);
		sin54sq = getSinSqTerm(minModelXSec.dm54Sq,E1,E2,l1,l2);
	}
	else{
		sin51sq = 0;
		sin54sq = 0;
	}

	if(minModelXSec.dm61Sq != 0){
		sin61sq = getSinSqTerm(minModelXSec.dm61Sq,E1,E2,l1,l2);
		sin64sq = getSinSqTerm(minModelXSec.dm64Sq,E1,E2,l1,l2);
		sin65sq = getSinSqTerm(minModelXSec.dm65Sq,E1,E2,l1,l2);
	}
	else{
		sin61sq = 0;
		sin64sq = 0;
		sin65sq = 0;
	}

	double osc_int = 4 * ((1 - pow(minModelXSec.Ue[0],2) - pow(minModelXSec.Ue[1],2) - pow(minModelXSec.Ue[2],2))
			* (pow(minModelXSec.Ue[0],2) * sin41sq + pow(minModelXSec.Ue[1],2) * sin51sq + pow(minModelXSec.Ue[2],2)*sin61sq)
			+ pow(minModelXSec.Ue[0],2) * pow(minModelXSec.Ue[1],2) * sin54sq
			+ pow(minModelXSec.Ue[0],2) * pow(minModelXSec.Ue[2],2) * sin64sq
			+ pow(minModelXSec.Ue[1],2) * pow(minModelXSec.Ue[2],2) * sin65sq);
  */

	return osc_int;
}

void fcnXSec(int &npar, double *gin, double &fval, double *xval, int iflag){

	double *_gin = gin;	int &_npar = npar;	int _iflag = iflag;

	int nBins_lsnd = 5;
	int nBins_karmen = 6;
	double d_karmen = 16.2;	double l_karmen = 3.;
	double d_lsnd = 25.475;	double l_lsnd = 8.75;

	double k_lsnd  = xval[0];
	double k_karmen = xval[1];
	double k_correl = xval[2];

	double chisq;
	chisq = pow((k_lsnd - 1.)/Lsnd_sys,2) + pow((k_karmen - 1.)/Karmen_sys,2) + pow((k_correl - 1.)/Correl_sys,2);

	double delE, lLo, lHi, ELo, EHi, EAvg, oscProb, xSec;
	// First, let's loop through karmen
	for(int i = 0; i < nBins_karmen; i++){
		if(i > 0)
			delE = (Karmen_Enu[i] - Karmen_Enu[i-1])/2.;
		else
			delE = (Karmen_Enu[i+1] - Karmen_Enu[i])/2.;
		lLo = d_karmen;
		lHi = d_karmen + l_karmen;
		ELo = Karmen_Enu[i] - delE;
		EHi = Karmen_Enu[i] + delE;
		EAvg = (EHi + ELo)/2.;
		ESigma = max(.08, .115/sqrt(EAvg - 17.3));
		oscProb = osc_int(ELo,EHi,lLo,lHi);
		xSec = (-2.5954e-4 * pow(EAvg,3) + 5.0028e-2 * pow(EAvg,2) - 1.5280 * EAvg + 12.876) * (1 - oscProb) * k_karmen * k_correl;
		chisq += pow((Karmen[i] - xSec)/Karmen_error[i],2);
	}

	// Now, do whatever we just did, but to LSND
	for(int i = 0; i < nBins_lsnd; i++){
		if(i > 0)
			delE = (Lsnd_Enu[i] - Lsnd_Enu[i-1])/2.;
		else
			delE = (Lsnd_Enu[i+1] - Lsnd_Enu[i])/2.;
		lLo = d_lsnd;
		lHi = d_lsnd + l_lsnd;
		ELo = Lsnd_Enu[i] - delE;
		EHi = Lsnd_Enu[i] + delE;
		EAvg = (EHi + ELo)/2.;
		ESigma = max(.08, .48/sqrt(EAvg - 17.3));
		oscProb = osc_int(ELo,EHi,lLo,lHi);
		xSec = (-2.5954e-4 * pow(EAvg,3) + 5.0028e-2 * pow(EAvg,2) - 1.5280 * EAvg + 12.876) * (1 - oscProb) * k_lsnd * k_correl;
		chisq += pow((Lsnd[i] - xSec)/Lsnd_error[i],2);
	}

	fval = chisq;

	return;
}
