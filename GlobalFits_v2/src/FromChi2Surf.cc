#include "FromChi2Surf.h"

int FromChi2Surf::Init(std::string dataLoc, Oscillator osc, bool debug){

  TFile * f = new TFile((dataLoc+"minosplus_dataRelease.root").c_str(),"READ");
  TH2D * surf = (TH2D*)f->Get("dm241vsth24");

  mySurf = (TH2D*)surf->Clone();
  mySurf->SetDirectory(0);

  f->Close();

  dof = 140;

  //Initialize output tree
  chi2Nt = new OutTree("MINOSplus");

  if(debug) std::cout << "MINOS+ initialized. Bins: " << 140 << std::endl;

  return dof;
}

float FromChi2Surf::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;

  double dm2, sin22th, theta;
  dm2 = pow(model.mNu[0],2);
  sin22th = 4 * pow(model.Ue[0],2) * pow(model.Um[0],2);
  theta = asin(sqrt(sin22th));

  if(theta > .01){
  //  //std::cout << "BOOOOM " << sin22th << " " << theta << " " << dm2 << std::endl;
    chi2 = mySurf->Interpolate(theta,dm2);
  }
  else chi2 = 99999;

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "Minos+ Chi2: " << chi2 << std::endl;
  return chi2;
}
