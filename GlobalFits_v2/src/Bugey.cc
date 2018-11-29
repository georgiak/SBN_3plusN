#include "Bugey.h"

// Minuit requires a static function to minimize, but I still need to give it inputs
// So i'm going to make some global variables here for it to take from. this isn't ideal, but it'll do.
double SigmaBigA = 4.796e-2;
double SigmaB =2.e-2;
double SigmaSmallA =1.414e-2;
bool ReactorAnomaly = false;
TMinuit *gMinuit;

double* dm2VecBugey;
std::vector < std::vector < double > > ObservedBugey, SigmaRatioBugey, EnergyBugey;
std::vector < std::vector < std::vector < double > > > sinSqDeltaGridBugey;
neutrinoModel minModelBugey;


int Bugey::Init(std::string dataLoc, Oscillator osc, bool debug){

  ObservedBugey.resize(nBaselines, std::vector<double>(maxEnergyBins));
  SigmaRatioBugey.resize(nBaselines, std::vector<double>(maxEnergyBins));
  EnergyBugey.resize(nBaselines, std::vector<double>(maxEnergyBins));
  sinSqDeltaGridBugey.resize(dm2VecMaxDim, std::vector<std::vector<double>>(maxEnergyBins, std::vector<double>(nBaselines)));

  dm2VecBugey = osc.dm2Vec;

  double deltaE[] = {.2, .2, .5};
  double EMin = 1.; double EMax = 6.;

  double ratio_obs1[] = {1.0182,1.0017,0.9815,0.9875,0.9981,0.9849,0.9749,0.9743,0.9840,0.9910,1.0042,0.9966,0.9629,1.0144,
      0.9971,0.9759,0.9889,0.9493,1.0056,0.8928,1.0160,0.9291,0.8585,1.0012,0.8981};
  double ratio_obs2[] = {0.9589,1.0344,0.9899,1.0071,0.9990,1.0267,0.9830,1.0298,0.9889,0.9441,0.9827,1.0123,1.0069,0.9059,
      0.9877,1.0281,1.0104,0.8905,1.0377,0.9940,1.0593,0.9457,0.9146,0.9941,0.9161};
  double ratio_obs3[] = {0.1685,0.8189,1.2169,1.3342,0.7995,1.1623,1.2918,1.3251,1.2711,0.5812};

  double sigmaRatio1[] = {0.01893,0.01610,0.01516,0.01472,0.01461,0.01463,0.01522,0.01595,0.01660,0.01906,0.02061,0.02064,
      0.02071,0.02366,0.02498,0.02681,0.02914,0.03117,0.03657,0.03699,0.04789,0.05124,0.05401,0.07408,0.07716};
  double sigmaRatio2[] = {0.06553,0.04876,0.03628,0.03182,0.03129,0.03140,0.03128,0.03082,0.02922,0.02847,0.03159,0.03412,
      0.03705,0.03669,0.04426,0.04753,0.05371,0.05415,0.06137,0.07084,0.07786,0.08788,0.10928,0.13213,0.17449};
  double sigmaRatio3[] = {0.54296,0.29210,0.20995,0.18166,0.23880,0.26919,0.34341,0.28928,0.79799,0.90417};

  double normReactorAno[] = {1.06237,1.06197,1.0627};
  double nBins[] = {25, 25, 10};

  double baselines[] = {15., 40., 95.};

  for(int j = 0; j < nBaselines; j++){
    for(int i = 0; i < nBins[j]; i++){
      EnergyBugey[j][i] = EMin + (double(i) + 0.5)/(nBins[j])*(EMax - EMin);
      if(j==0){
        ObservedBugey[j][i] = ratio_obs1[i];
        SigmaRatioBugey[j][i] = sigmaRatio1[i];
      }
      if(j==1){
        ObservedBugey[j][i] = ratio_obs2[i];
        SigmaRatioBugey[j][i] = sigmaRatio2[i];
      }
      if(j==2){
        ObservedBugey[j][i] = ratio_obs3[i];
        SigmaRatioBugey[j][i] = sigmaRatio3[i];
      }
    }
  }

  // Now, for the more complicated stuff
  double EnuMin = 1.8; double EnuMax = 10.;
  integralFuncsBugey integralFuncs;
  ROOT::Math::Functor1D wf(&integralFuncs, &integralFuncsBugey::sinSqFunction);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ig.SetAbsTolerance(1*pow(10,-4));
  ig.SetRelTolerance(1*pow(10,-4));
  ig.SetFunction(wf);

  for(int k = 0; k < dm2VecMaxDim; k++){
    for(int j = 0; j < nBaselines; j++){
      for(int i = 0; i < nBins[j]; i++){

        double Enu = EnergyBugey[j][i] + 1.8;
        int n = 1.27 * dm2VecBugey[k] * baselines[j] / (Enu * TMath::Pi());
        double eNuNext = 1.27 * dm2VecBugey[k] * baselines[j] / ((n + 0.5) * TMath::Pi());

        if(abs(Enu - eNuNext) >= 0.123 * sqrt(Enu - 1.8)){
          integralFuncs._energy = Enu;
          integralFuncs._dm2 = dm2VecBugey[k];
          integralFuncs._jB = j;

          sinSqDeltaGridBugey[k][i][j] = ig.Integral(EnuMin,EnuMax);// / (EnuMax-EnuMin);
        }
        else
          sinSqDeltaGridBugey[k][i][j] = .5;
      }
    }
  }

  // Initialize Minuit for fitting later
  gMinuit = new TMinuit(6);

  dof = 60;

  //Initialize output tree
  chi2Nt = new  OutTree("Bugey");

  if(debug) std::cout << "Bugey initialized. Bins: " << dof << std::endl;
  return dof;
}

float Bugey::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;

  minModelBugey = model;
  double chisq = 0;

  double arglis[2];
  int ierflag;
  arglis[0] = -1.;
  gMinuit->SetFCN(fcnBugey);
  gMinuit->mnexcm("SET PRI",arglis,1,ierflag);
  // Okay, let's get this minuit garbage started
  arglis[0] = 1.;
  gMinuit->mnexcm("SET STR",arglis,1,ierflag);
  // Now, clear parameters and set 'em anew
  gMinuit->mnexcm("CLE",arglis,1,ierflag);
  gMinuit->mnparm(0,TString("bigA"),1.,SigmaBigA,.5,1.5,ierflag);
  gMinuit->mnparm(1,TString("smallA1"),1.,SigmaSmallA,.5,1.5,ierflag);
  gMinuit->mnparm(2,TString("smallA2"),1.,SigmaSmallA,.5,1.5,ierflag);
  gMinuit->mnparm(3,TString("smallA3"),1.,SigmaSmallA,.5,1.5,ierflag);
  gMinuit->mnparm(4,TString("b"),0.,SigmaB,-.1,.1,ierflag);

  // Now reset the function call and check if we even need to minimize. If it's a lost cause, don't waste the computation time.
  arglis[0] = 5.;
  gMinuit->mnexcm("CALL",arglis,1,ierflag);
  double disttomin, errdef;   int npari, nparx, istat;

  arglis[0] = 10000.;
  arglis[1] = 100.;
  gMinuit->mnexcm("MINI",arglis,2,ierflag);
  gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
  chi2 = chisq;

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "Bugey Chi2: " << chi2 << std::endl;
  return chi2;
}

double  integralFuncsBugey::sinSqFunction(const double x){

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

void fcnBugey(int &npar, double *gin, double &fval, double  *xval, int iflag){
  // This function will calculate the chi2. it's to be minimized with MINUIT,  which is why it's so convoluted.
	double *_gin = gin;	int &_npar = npar;	int _iflag = iflag;

  double normReactorAno[] = {1.06237,1.06197,1.0627};
  double nBins[] = {25, 25, 10};
  double bigA, b, smallA[3];

  bigA = xval[0];
  smallA[0] = xval[1];
  smallA[1] = xval[2];
  smallA[2] = xval[3];
  b = xval[4];

  double prob, chisq;
  double Ei = 1.;
  double sinSq[dm2VecMaxDim];
  ROOT::Math::Interpolator dif(dm2VecMaxDim,ROOT::Math::Interpolation::kCSPLINE);
	oscContribution oscCon = getOscContributionsNueDis(minModelBugey);

  chisq = 0.;

  for(int j = 0; j < 3; j++){ // loop over baselines
	  for(int i = 0; i < nBins[j]; i++){  // loop over energy bins
  	  for(int k = 0; k < dm2VecMaxDim; k++){
        sinSq[k] = sinSqDeltaGridBugey[k][i][j];
      }
      // We're doing nue disappearance
			prob = 1.;
      dif.SetData(dm2VecMaxDim,dm2VecBugey,sinSq);
      for(int iCon = 0; iCon < 6; iCon ++){
			if(oscCon.dm2[iCon] != 0.){
					prob += oscCon.aEE[iCon] * dif.Eval(oscCon.dm2[iCon]);
				}
			}

      // Now, calculate the chisq
      double num = (bigA * smallA[j] + b * (EnergyBugey[j][i] - Ei)) * prob - ObservedBugey[j][i];
			if(ReactorAnomaly)
    			num = (bigA * smallA[j] + b * (EnergyBugey[j][i] - Ei)) * prob * normReactorAno[j] - ObservedBugey[j][i];
			chisq += pow(num/SigmaRatioBugey[j][i],2);
    	}
    	chisq += pow((smallA[j] - 1.)/SigmaSmallA,2);
  	}

  	chisq += pow((bigA - 1.)/SigmaBigA,2) + pow(b/SigmaB,2);
  	fval = chisq;
}
