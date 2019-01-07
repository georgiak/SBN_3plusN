#include "PROSPECT.h"

/*
Davio Cianci Oct 15, 2018

PROSPECT paper: https://arxiv.org/abs/1806.02784

Current fit just uses chi2 surface. Note, we ignore critical chi2s given in paper and do non-frequentist way
assuming no correlation

*/

int PROSPECT::Init(std::string dataLoc, Oscillator osc, bool debug){

  surf = new TGraph2D();

  ifstream file;
  // Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
  file.open(dataLoc+"prospect_chi2surf.txt");
  double sin22th, dm2, deltachi2;
  for(int i = 0; i < 2520; i++){
    file >> sin22th;
    file >> dm2;
    file >> deltachi2;
    surf->SetPoint(i,sin22th,dm2,deltachi2);
  }

  file.close();

  dof = 80;

  // Iniitalize our output tree
  chi2Nt = new OutTree("PROSPECT");

	if(debug) std::cout << "PROSPECT initialized. Bins: " << 80 << std::endl;

  return dof;
}

float PROSPECT::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

	oscContribution oscCon = getOscContributionsNueDis(model);

	double dm2, sin22th;
  double chi2min = 57.9;
  dm2 = pow(model.mNu[0],2);
	sin22th = 4 * pow(model.Ue[0],2) * (1 - pow(model.Ue[0],2));
  chi2 = surf->Interpolate(sin22th,dm2) + chi2min;

  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "PROSPECT Chi2: " << chi2 << std::endl;

  return chi2;
}
