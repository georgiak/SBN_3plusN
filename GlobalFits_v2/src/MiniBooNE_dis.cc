#include "MiniBooNE_dis.h"

int MiniBooNE_dis::Init(std::string dataLoc, Oscillator osc, bool debug){

  std::string str_data_numu, str_fracterrormatrix, str_fullosc, str_binboundaries;
  str_binboundaries = dataLoc + "miniboone_dis/miniboone_binboundaries_disap.txt";

  Signal.resize(nBins);
  FullData.resize(nBins);
  Full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));

  if(!nubar){
    // If in neutrino mode
    nFOsc = 1267007;
    str_data_numu = dataLoc + "miniboone_dis/miniboone_numudata_disap.txt";
    str_fracterrormatrix = dataLoc + "miniboone_dis/miniboone_frac_shape_matrix_numu_disap.txt";
    str_fullosc = dataLoc + "miniboone_dis/numudisap_ntuple.txt";
  }
  else{
    // If in antineutrino mode
    nFOsc = 686529;
    str_data_numu = dataLoc + "miniboone_dis/miniboone_numubardata_disap.txt";
    str_fracterrormatrix = dataLoc + "miniboone_dis/miniboone_frac_shape_matrix_numubar_disap.txt";
    str_fullosc = dataLoc + "miniboone_dis/numubardisap_ntuple.txt";
  }

  float *nu_EnuQE = new float[nBins + 1];
	float *nu_FOsc_EnuQE = new float[nFOsc];
	float *nu_FOsc_EnuTrue = new float[nFOsc];
	float *nu_FOsc_LnuTrue = new float[nFOsc];
	float *nu_FOsc_weight = new float[nFOsc];

  ifstream file;
  file.open(str_data_numu);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << str_data_numu << std::endl;
	for(short i = 0; i < nBins; i++)
		file >> FullData[i];
	file.close();

	file.open(str_fracterrormatrix);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << str_fracterrormatrix << std::endl;
	for(short i = 0; i < nBins; i++)
		for(short j = 0; j < nBins; j++){
			file >> Full_fractCovMatrix[i][j];
    }
	file.close();

  file.open(str_binboundaries);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << str_binboundaries << std::endl;
  for(short i = 0; i < nBins+1; i++)
      file >> nu_EnuQE[i];
  file.close();

 	file.open(str_fullosc);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << str_fullosc << std::endl;
	int dummy;
    for(int iEvt = 0; iEvt < nFOsc; iEvt++){
      file >> dummy;
      file >> nu_FOsc_EnuQE[iEvt];
      file >> nu_FOsc_EnuTrue[iEvt];   // true energy of neutrino
      file >> nu_FOsc_LnuTrue[iEvt];   // distance from production and detection points
      file >> nu_FOsc_weight[iEvt];    // event weight
	}
	file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
	float mstep = TMath::Log10(100./.01)/float(100);
	Libdis_sinsq.resize(100, std::vector<float>(nBins));
	Libdis_noosc.resize(nBins);

	for(int mi = 0; mi < 100; mi++){
		float dm2 = pow(10,((mi+1.)/100.*TMath::Log10(100./.01) + TMath::Log10(.01)));
		for(int iB = 0; iB < nBins; iB++){
			Libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				Libdis_noosc[iB] = 0;
		}

		for(int iFOsc = 0; iFOsc < nFOsc; iFOsc++){   // Loop over full oscillation events
      for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	      if(nu_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nu_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){

					float ETru = nu_FOsc_EnuTrue[iFOsc];
					float LTru = nu_FOsc_LnuTrue[iFOsc];
					Libdis_sinsq[mi][iB] += nu_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru/ETru),2); // oscillated events in each bin
					if(mi == 0)
						Libdis_noosc[iB] += nu_FOsc_weight[iFOsc];  // unoscillated events in each bin
        }
      }
    }
  }
  dof = nBins;

  // Initialize output tree
  if(!nubar)
    chi2Nt = new OutTree("MBnu_dis");
  else
    chi2Nt = new OutTree("MBnubar_dis");

  if(debug){
    if(!nubar) std::cout << "MBnu_dis initialized. Bins: " << nBins << std::endl;
    else std::cout << "MBnubar_dis initialized. Bins: " << nBins << std::endl;
  }
  return dof;

}

float MiniBooNE_dis::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;
  std::vector <  float > Prediction, PredictionNoOsc;
  Prediction.resize(nBins);
  PredictionNoOsc.resize(nBins);

  // Initialize contributions from osc probability
	double sin22th = model.ProbAmp("mumu");

  covMatrix.ResizeTo(nBins, nBins);
  covMatrix.Zero();

  float minEBins[nBins], maxEBins[nBins];
  float dtIntegral(0.), MCIntegral(0.), MCIntegralNoOsc(0.);

	for(int iB = 0; iB < nBins; iB ++){
		Signal[iB] = 0;
		PredictionNoOsc[iB] = Libdis_noosc[iB];
		dtIntegral += FullData[iB];
    MCIntegralNoOsc += Libdis_noosc[iB];
	}

	float mstep = TMath::Log10(100./.01)/float(100);
	int dm2;
	for(int iB = 0; iB < nBins; iB++){
		dm2 = floor(TMath::Log10(model.Dm2()/.01)/mstep);
    Prediction[iB] = PredictionNoOsc[iB] - sin22th*Libdis_sinsq[dm2][iB];
	}
	// Normalize signal prediction to data
	for(int iB = 0; iB < nBins; iB++){
		MCIntegral += (Prediction[iB] + Signal[iB]);
	}

	for(int iB = 0; iB < nBins; iB++){
		Prediction[iB] *= dtIntegral/MCIntegral;
    PredictionNoOsc[iB] *= dtIntegral/MCIntegralNoOsc;
	}

  //std::cout << "Integrals: " << dtIntegral << " " << MCIntegral << std::endl;

	for(int iB = 0; iB < nBins; iB++){
		for(int jB = 0; jB < nBins; jB++){
			covMatrix[iB][jB] = Full_fractCovMatrix[iB][jB] * Prediction[iB] * Prediction[jB];
			// Add statistical error of signal prediction
			if(iB == jB){
				covMatrix[iB][jB] += Prediction[iB];
        //std::cout << "stat err: " << Prediction[iB] << std::endl;
			}
		}
	}
	// Now, let's invert this bad boy
	cov.ResizeTo(nBins,nBins);
	cov = covMatrix.Invert();

	// Finally, let's put everything together and calculate the chisq
	for(int iB = 0; iB < nBins; iB++){
		for(int jB = 0; jB < nBins; jB++){
			chi2 += (FullData[iB] - Prediction[iB]) * cov[iB][jB] * (FullData[jB] - Prediction[jB]);
		}
	}

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    if(!nubar)  std::cout << "MBnu_dis Chi2: " << chi2 << std::endl;
    else std::cout << "MBnubar_dis Chi2: " << chi2 << std::endl;
  }
  return chi2;
}
