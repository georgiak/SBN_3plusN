#include "IceCube_fromSurfs.h"

int IceCube_fromSurfs::Init(std::string dataLoc, Oscillator osc, bool debug){

  //
  // Ok. Icecube is an ordeal, so I put this together mainly for validation studies.
  // All of the icecube work is done elsewhere and packed into a nice array of 3d histograms.
  // Here, we read those histograms and properly interpolate them.
  //

  TFile * f = new TFile((dataLoc+"IC_chi2grids.root").c_str(),"READ");
  for(int i = 0; i < 10; i++){
    vecLogICHisto.push_back((TH3D*)f->Get(("loghisto_"+std::to_string(i)).c_str())->Clone());
    vecLogICHisto[i]->SetDirectory(0);
  }
  f->Close();

  dof = 300;

  //Initialize output tree
  chi2Nt = new OutTree("IceCube");

  if(debug) std::cout << "IceCube initialized. Bins: " << dof << std::endl;

  return dof;
}

float IceCube_fromSurfs::Chi2(Oscillator osc, neutrinoModel model, bool debug){
  double chi2 = 99999;

  // Ok. So we have an array of 3d histograms.
  // Each TH3D has a chi2 surface slice in theta34 with dimensions log10(m41)x log10(th14)xlog10(th24)
  // Logs are so we can do a smooth linear interpolation
  // We want to find the MIN, which means looping through all 10 histograms at the same point

  std::array<double,4> ops = model.OscParams();

  if(ops[0] < .1 || ops[0] > 10.0 || ops[1] < .01 || ops[1] > 3.14/4 || ops[2] < .01 || ops[2] > 3.14/4){
    std::cout << "Osc Params out of bounds for icecube interpolation. Returning big number" <<  std::endl;
    chi2 = 99999;
  }
  else{
    for(int i = 0; i < 10; i++){
      chi2 = std::min(vecLogICHisto[i]->Interpolate(log10(ops[0]),log10(ops[1]),log10(ops[2])),chi2);
    }
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "IC Chi2: " << chi2 << std::endl;
  return chi2;
}
