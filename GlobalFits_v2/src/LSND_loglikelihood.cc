#include "LSND_loglikelihood.h"

int LSND_loglikelihood::Init(std::string dataLoc, Oscillator osc, bool debug){

  dm2Vec.resize(dm2VecMaxDim);
  sinSqDeltaGrid.resize(dm2VecMaxDim, std::vector<double>(nBins));
  sinSqDeltaGrid2.resize(dm2VecMaxDim, std::vector<double>(nBins));
  Bkg.resize(nBins);
  Observed.resize(nBins);

  // Load up data
  double temp0[] = {10., 11.7, 17.6, 7.8, 3.7};
  double temp1[] = {5.1, 5.2, 3.7, 2.2, 0.6};
  for(int iL = 0; iL < nBins; iL++){
    Observed[iL] = temp0[iL];
    Bkg[iL] = temp1[iL];
  }

  norm = 57401.252; // Taken from old code since multiple integrals in root are tricky and time consuming

  ifstream file;
  file.open(dataLoc+"lsndsinsq.txt");
	for(int k = 0; k < 601; k++){
        for(int iL = 0; iL < nBins; iL++){
			file >> dm2Vec[k];
			file >> sinSqDeltaGrid[k][iL];
			file >> sinSqDeltaGrid2[k][iL];
		}
	}
	file.close();

  dof = nBins;

  //Initialize output tree
  chi2Nt = new OutTree("LSND_loglikelihood");

  if(debug) std::cout << "LSND initialized. Bins: " << nBins << std::endl;

  return dof;
}

float LSND_loglikelihood::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;

  double sin22th = model.ProbAmp("mue");
  double dm2 = model.Dm2();

  // Initialize Inerpolator
  ROOT::Math::Interpolator dif(dm2VecMaxDim);
  double sinSq, sinSq2;
  double sinSqDeltaVec[dm2VecMaxDim], sinSqDeltaVec2[dm2VecMaxDim], lsndDm2Vec[dm2VecMaxDim];

  std::vector < double > Signal;
  Signal.assign(nBins,0.);

  for(int iL = 0; iL < nBins; iL ++){
    for(int k = 0; k < 601; k ++){
      sinSqDeltaVec[k] = sinSqDeltaGrid[k][iL];
      sinSqDeltaVec2[k] = sinSqDeltaGrid2[k][iL];
      lsndDm2Vec[k] = dm2Vec[k];
    }
    dif.SetData(dm2VecMaxDim,lsndDm2Vec,sinSqDeltaVec);
    Signal[iL] += sin22th * dif.Eval(dm2) * norm;
  }

  // Now, using the signal vector, use the log-likelihood method to get the effective chisq
  double lt1, lt2, lt3, pred;
  for(int iL = 0; iL < nBins; iL++){
      pred = Signal[iL] + Bkg[iL];
      lt1 = pred - Observed[iL];
      if(pred > 0) lt2 = Observed[iL] * log(pred);
      else lt2 = 0;
      if(Observed[iL] > 0) lt3 = Observed[iL] * log(Observed[iL]);
      else lt3 = 0;

      chi2+= 2 * (lt1 - lt2 + lt3);
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "LSND Chi2: " << chi2 << std::endl;
  return chi2;
}
