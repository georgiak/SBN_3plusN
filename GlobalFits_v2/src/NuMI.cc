#include "NuMI.h"

int NuMI::Init(std::string dataLoc, Oscillator osc, bool debug){

  Signal.resize(nBins);
  TotalError.resize(nBins);
  NueData.resize(nBins);
  NueBgr.resize(nBins);
  NueBgr_error.resize(nBins);
  FOsc_fracError.resize(nBins);

  float *EnuQE = new float[nBins+1];
  float *FOsc_EnuQE = new float[nFOscEvts];     // reconstructed neutrino energy
  float *FOsc_EnuTrue = new float[nFOscEvts];   // true energy of neutrino
  float *FOsc_LnuTrue = new float[nFOscEvts];   // distance from production and detection points
  float *FOsc_weight = new float[nFOscEvts];

  float temp0[] = {200.,300.,475.,675.,900.,1100.,1300.,1500.,1700.,2000.,3000.};
  int temp1[] = {59,142,151,146,83,68,57,39,19,16};
  float temp2[] = {41.7574,117.537,123.673,118.025,82.0508,64.8166,41.1826,31.3267,22.0301,17.6672};
  float temp3[] = {12.112,25.5363,27.4975,24.8745,18.1301,15.11,10.5605,9.03036,6.88422,5.70684};

  for(int i = 0; i < nBins; i++){
    EnuQE[i] = temp0[i];
    NueData[i] = temp1[i];
    NueBgr[i] = temp2[i];
    NueBgr_error[i] = temp3[i];
  }
  EnuQE[nBins] = temp0[nBins];

  ifstream file;
  // Get nFullOscEvts from another file!
  file.open(dataLoc+"numi_fullosc.out");
  for(int i = 0; i < nFOscEvts; i++){
      file >> FOsc_weight[i];     // reconstructed neutrino energy
      file >> FOsc_EnuTrue[i];   // true energy of neutrino
      file >> FOsc_EnuQE[i];   // distance from production and detection points
      file >> FOsc_LnuTrue[i];    // event weight
  }
  file.close();

  for(int i = 0; i < nBins; i++){
      FOsc_fracError[i] = NueBgr_error[i] / NueBgr[i];
  }

  // Precompute sines and sinesquareds
  float dmmax = 100.;
  float dmmin = 0.01;
  float mstep = TMath::Log10(dmmax/dmmin)/100.f;
  float ETru, LTru, dm2;
  Lib_sinsq.resize(100, std::vector<float>(nBins));
  Lib_sin.resize(100, std::vector <float>(nBins));

  for(int mi = 0; mi < 100; mi++){

    dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
    for(int iB = 0; iB < nBins; iB++){
      Lib_sinsq[mi][iB] = 0;
      Lib_sin[mi][iB] = 0;
    }
    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
      for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

        if(FOsc_EnuQE[iFOsc] > EnuQE[iB] && FOsc_EnuQE[iFOsc] < EnuQE[iB+1]){
          ETru = FOsc_EnuTrue[iFOsc];
          LTru = FOsc_LnuTrue[iFOsc];

          Lib_sinsq[mi][iB] += FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
          Lib_sin[mi][iB] += FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
        }
      }
    }
  }

  dof = nBins;

  //Initialize output tree
  chi2Nt = new OutTree("NuMI");

  if(debug) std::cout << "NuMI initialized. Bins: " << nBins << std::endl;

  return dof;
}

float NuMI::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;
  std::vector < float > Prediction;
  Prediction.resize(nBins);

  // Initialize contributions from the osc probability
  oscContribution oscCont;
  oscCont = getOscContributionsNueApp(model, false, true);

  Signal.assign(nBins,0.);
  TotalError.assign(nBins,0.);

  float mstep = TMath::Log10(100./.01)/float(100);
  int dm2;
  for(int iB = 0; iB < nBins; iB++){
  	for(int iContribution = 0; iContribution < 6; iContribution++){
  		if(oscCont.dm2[iContribution] == 0)	Signal[iB] += 0;
  		else{
  			dm2 = floor(TMath::Log10(oscCont.dm2[iContribution]/.01)/mstep);
  			Signal[iB] += oscCont.aMuE[iContribution]*Lib_sinsq[dm2][iB] + oscCont.aMuE_CPV[iContribution]*Lib_sin[dm2][iB];
  		}
  	}
  }
  // Now fill up that error
  for(int iN = 0; iN < nBins; iN++){
    TotalError[iN] = pow(NueBgr_error[iN],2) + pow(FOsc_fracError[iN]*Signal[iN],2) + Signal[iN] + NueBgr[iN];
  }

  // Now actually calculate the chisq
  for(int iN = 0; iN < nBins; iN++){
    chi2 += pow((NueData[iN] - NueBgr[iN] - Signal[iN]),2) / TotalError[iN];
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "NuMI Chi2: " << chi2 << std::endl;
  return chi2;
}
