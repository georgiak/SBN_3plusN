#include "Atm.h"

int Atm::Init(std::string dataLoc, Oscillator osc, bool debug){

  DmuVec = new double[dmuVecMaxDim];
  DChi2Vec = new double[dmuVecMaxDim];
  double dmuMin = 0.;
  double dmuMax = .25;

  double temp0[] = {0.00000,0.03913,0.10552,0.18864,0.28771,0.40839,0.54925,0.69901,0.86038,1.05169,1.27231,1.53027,1.77613,2.02014,2.27515,2.54545,2.83366,3.13347,3.44052,3.76300,
              4.09784,4.46017,4.84822,5.26111,5.69840,6.16382,6.65471,7.16706,7.70432,8.24987,8.80458,9.38657,9.99066,10.61653,11.26452,11.93572,12.62999,13.32833,14.03633,14.77024,
              15.52484,16.29127,17.07689,17.87941,18.69526,19.53419,20.39532,21.27474,22.17658,23.08987,24.01778,24.969511,25.93752,26.92796,27.94275,28.98161,30.04210,31.12363,32.22914,33.35148,
              34.49207,35.65003,36.82462,38.01945,39.23037,40.46261,41.71198,42.97714,44.26171,45.56482,46.88852,48.22036,49.56687,50.93319,52.30869,53.58208,54.86813,56.16923,57.48430,58.81740,
              60.16430,61.52579,62.90204,64.28466,65.67754,67.07690,68.48444,69.90443,71.33371,72.77951,74.23951,75.70818,77.18138,78.66365,80.15671,81.66008,83.17423,84.69586,86.22622,87.77160,89.32895};

  for(int iA = 0; iA < dmuVecMaxDim; iA++){
      DmuVec[iA] = dmuMin + (dmuMax - dmuMin) * (iA)/float(dmuVecMaxDim-1);
      DChi2Vec[iA] = temp0[iA];
  }

  dof = 1;

  //Initialize output tree
  chi2Nt = new OutTree("Atm");

  if(debug) std::cout << "Atm initialized. Bins: " << 1 << std::endl;

  return dof;
}

float Atm::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;

  ROOT::Math::Interpolator dif(dmuVecMaxDim);
  dif.SetData(dmuVecMaxDim,DmuVec,DChi2Vec);

  // First, let's find dmu such that dmu**2 - dmu + A = 0
  //double A = (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2))  * (pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2)) + pow(model.Um[0]*model.Um[1],2) + pow(model.Um[0]*model.Um[2],2) + pow(model.Um[1]*model.Um[2],2);
  //double dmu = (1 - sqrt(1. - 4*A)) /2.;

  // Ok. According to https://arxiv.org/pdf/hep-ph/0405172.pdf, the following formula for dmu applies for 3+1 but will need to be generalized for 3+N
  double dmu = 1 - pow(model.Umu4(),2);
  double A = (1. - pow(model.Umu4(),2))  * (pow(model.Umu4(),2));
  double dmutwo = (1 - sqrt(1. - 4*A)) /2.;


  std::cout << "dmu = " << dmu << " " << dmutwo << std::endl;

  // For this, everything has been figured out, so we just interpolate an array of chi2's to get our result as a function of dmu
  chi2 = dif.Eval(dmu);

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "Atm Chi2: " << chi2 << std::endl;
  return chi2;
}
