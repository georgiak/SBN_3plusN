   #include "MiniBooNE.h"

int MiniBooNE::Init(std::string dataLoc, bool debug){

  std::string str_data_nue, str_data_numu, str_MC_nue, str_MC_numu, str_fracterrormatrix, str_binboundaries, str_fullosc;
  str_binboundaries = dataLoc + "miniboone_binboundaries_lowe.txt";

  Signal.resize(nBins_e);
  FullData.resize(nBins_e + nBins_mu);
  Background.resize(nBins_e + nBins_mu);
  Full_fractCovMatrix.resize(nBins_e + nBins_e + nBins_mu, std::vector < float > (nBins_e + nBins_e + nBins_mu));

  float *nu_EnuQE = new float[nBins_e + 1];
	float *nu_FOsc_EnuQE = new float[nFOscEvts];
	float *nu_FOsc_EnuTrue = new float[nFOscEvts];
	float *nu_FOsc_LnuTrue = new float[nFOscEvts];
	float *nu_FOsc_weight = new float[nFOscEvts];

  if(!nubar){
    // If in neutrino mode:
    str_data_nue = dataLoc + "miniboone/data_nue.txt";
    str_data_numu = dataLoc + "miniboone/data_numu.txt";
    str_MC_nue = dataLoc + "miniboone/MC_nue.txt";
    str_MC_numu = dataLoc + "miniboone/MC_numu.txt";
    str_fracterrormatrix = dataLoc + "miniboone/frac_error_matrix_contNu.txt";
    str_fullosc = dataLoc + "miniboone/numunuefullosc_ntuple.txt";
  }
  else{
    // If in antineutrino mode:
    str_data_nue = dataLoc + "miniboone_nuebardata_lowe.txt";
    str_data_numu = dataLoc + "miniboone_numubardata.txt";
    str_MC_nue = dataLoc + "miniboone_nuebarbgr_lowe.txt";
    str_MC_numu = dataLoc + "miniboone_numubar.txt";
    str_fracterrormatrix = dataLoc + "miniboone_full_fractcovmatrix_nubar_lowe";
    str_fullosc = dataLoc + "miniboone_numubarnuebarfullosc_ntuple.txt";
  }

  ifstream file;
  // Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
  file.open(str_data_nue);
  for(int i = 0; i < nBins_e; i++)
    file >> FullData[i];
  file.close();

  // Get measured Numu ccqe events per enuqe bin
  file.open(str_data_numu);
  for(int i = 0; i < nBins_mu; i++)
    file >> FullData[i + nBins_e];
  file.close();

  // Get predicted nue background events per enuqe bin
  file.open(str_MC_nue);
  for(int i = 0; i < nBins_e; i++)
    file >> Background[i];
  file.close();

  // Get predicted numu ccqe events per enuqe bin
  file.open(str_MC_numu);
  for(int i = 0; i < nBins_mu; i++)
    file >> Background[i + nBins_e];
  file.close();

  // Get fractional cov matrix for full numu->nue oscillation events
  file.open(dataLoc+"neutrino_frac_error_matrix.txt");
  for(int i = 0; i < nBins_e + nBins_e + nBins_mu; i++)
    for(int j = 0; j < nBins_e + nBins_e + nBins_mu; j++)
      file >> Full_fractCovMatrix[i][j];
  file.close();

  file.open(str_binboundaries);
  for(int i = 0; i < nBins_e + 1; i++)
    file >> nu_EnuQE[i];
  file.close();

  // Get nFullOscEvts for full numu->nue osc events after nue cuts
  file.open(str_fullosc);
  for(int iEvt = 0; iEvt < nFOscEvts; iEvt++){
        file >> nu_FOsc_EnuQE[iEvt];
        file >> nu_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nu_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nu_FOsc_weight[iEvt];    // event weight
  }
  file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
  float dmmax = 100.;
  float dmmin = 0.01;
  float mstep = TMath::Log10(dmmax/dmmin)/100.f;
  float ETru, LTru, dm2;
  Lib_sinsq.resize(100, std::vector<float>(nBins_e));
  Lib_sin.resize(100, std::vector<float>(nBins_e));

  for(int mi = 0; mi < 100; mi++){

    dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
    for(int iB = 0; iB < nBins_e; iB++){
      Lib_sinsq[mi][iB] = 0;
      Lib_sin[mi][iB] = 0;
    }
    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
      for(int iB = 0; iB < nBins_e; iB++){    // Loop over energy bins to fill the prediction vector pred

        if(nu_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nu_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){
          ETru = nu_FOsc_EnuTrue[iFOsc];
          LTru = nu_FOsc_LnuTrue[iFOsc];

          Lib_sinsq[mi][iB] += nu_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
          Lib_sin[mi][iB] += nu_FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
        }
      }
    }
  }

	// Lastly, Initialize best fit signal, with which we will weigh our covariance matrix
  neutrinoModel model_bestFit;
  model_bestFit.zero();
  model_bestFit.Ue[0] = 0.34;  model_bestFit.Um[0] = 0.34;    model_bestFit.mNu[0] = .19;
 	oscContribution oscCont_bestFit = getOscContributionsNueApp(model_bestFit, nubar, true);

	Signal_BestFit.resize(nBins_e + nBins_mu);

	for(int iB = 0; iB < nBins_e; iB++){
		for(int iContribution = 0; iContribution < 6; iContribution++){
    	if(oscCont_bestFit.dm2[iContribution] == 0) Signal_BestFit[iB] += 0;
			else{
      	dm2 = floor(TMath::Log10(oscCont_bestFit.dm2[iContribution]/.01)/mstep);
        Signal_BestFit[iB] += oscCont_bestFit.aMuE[iContribution]*Lib_sinsq[dm2][iB] + oscCont_bestFit.aMuE_CPV[iContribution]*Lib_sin[dm2][iB];
			}
		}
	}

  dof = nBins_e + nBins_mu - 1;

  // Initialize output tree
  if(!nubar)
    chi2Nt = new OutTree("MBnu");
  else
    chi2Nt = new OutTree("MBnubar");

  if(debug){
    if(!nubar) std::cout << "MBnu initialized. Bins: " << nBins_e + nBins_mu - 1 << std::endl;
    else std::cout << "MBnubar initialized. Bins: " << nBins_e + nBins_mu - 1 << std::endl;
  }
  return dof;
}

float MiniBooNE::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;
  std::vector < float > Prediction;
  Prediction.resize(nBins_e + nBins_mu);

  // Initialize contributions from the oscillation probability
  oscContribution oscCont;
  oscCont = getOscContributionsNueApp(model, nubar, true);

  // Zero out our signal
  std::fill(Signal.begin(), Signal.end(), 0);

  full_covMatrix.ResizeTo(nBins_e + nBins_e + nBins_mu, nBins_e + nBins_e + nBins_mu);
  full_covMatrix.Zero();
  covMatrix.ResizeTo(nBins_e + nBins_mu, nBins_e + nBins_mu);
  covMatrix.Zero();

  float mstep = TMath::Log10(100./.01)/float(100);
  int dm2;
  for(int iB = 0; iB < nBins_e; iB++){
  	for(int iContribution = 0; iContribution < 6; iContribution++){
  		if(oscCont.dm2[iContribution] == 0)	Signal[iB] += 0;
  		else{
  			dm2 = floor(TMath::Log10(oscCont.dm2[iContribution]/.01)/mstep);
  			Signal[iB] += oscCont.aMuE[iContribution]*Lib_sinsq[dm2][iB] + oscCont.aMuE_CPV[iContribution]*Lib_sin[dm2][iB];
  		}
  	}
  }
  // Divide signal prediction by the number of fullosc events
  for(int iB = 0; iB < nBins_e; iB++){
  	Signal[iB] /= float(nFOscEvts);
		Signal_BestFit[iB] /= float(nFOscEvts);
  }

  // Now, scale the fractional cov matrix to our signal and prediction vectors
  for(int iB = 0; iB < nBins_e + nBins_e + nBins_mu; iB++){
  	for(int jB = 0; jB < nBins_e + nBins_e + nBins_mu; jB++){
  		if(iB < nBins_e && jB < nBins_e){
  			full_covMatrix(iB,jB) = Full_fractCovMatrix[iB][jB]*Signal_BestFit[iB]*Signal_BestFit[jB];
  			// Add Stat error of signal prediction
  			if(iB == jB){
  				full_covMatrix(iB,jB) += Signal_BestFit[iB];
  			}
  		}
  		else if(iB < nBins_e && jB >= nBins_e){
  			full_covMatrix(iB,jB) = Full_fractCovMatrix[iB][jB]*Signal_BestFit[iB]*Background[jB-nBins_e];
  		}
  		else if(iB >= nBins_e && jB < nBins_e){
  			full_covMatrix(iB,jB)= Full_fractCovMatrix[iB][jB]*Background[iB-nBins_e]*Signal_BestFit[jB];
  		}
  		else if(iB >= nBins_e && jB >= nBins_e){
  			full_covMatrix(iB,jB) = Full_fractCovMatrix[iB][jB]*Background[iB-nBins_e]*Background[jB-nBins_e];
  		}
  	}
  }

  // Now, collapse our 3x3 matrix to a 2x2
  for(int iB = 0; iB < nBins_e + nBins_mu; iB++){
  	for(int jB = 0; jB < nBins_e + nBins_mu; jB++){
  		if(iB < nBins_e && jB < nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB][jB] + full_covMatrix[iB + nBins_e][jB] + full_covMatrix[iB][jB + nBins_e] + full_covMatrix[iB + nBins_e][jB + nBins_e];
  		}
  		else if(iB < nBins_e && jB >= nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB][jB + nBins_e] + full_covMatrix[iB + nBins_e][jB + nBins_e];
  		}
  		else if(iB >= nBins_e && jB < nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB + nBins_e][jB] + full_covMatrix[iB + nBins_e][jB + nBins_e];
  		}
  		else if(iB >= nBins_e && jB >= nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB + nBins_e][jB + nBins_e];
  	  }
  	}
  }

  // Now, let's invert the covariance matrix
  cov.ResizeTo(nBins_e + nBins_mu, nBins_e + nBins_mu);
  cov = covMatrix.Invert();

  for(int iB = 0; iB < nBins_e + nBins_mu; iB++)
    Prediction[iB] = Background[iB];

  for(int iB = 0; iB < nBins_e; iB++){
  	Prediction[iB] += Signal[iB];
  }

  // Finally, let's put everything together and calculate the chisq
  for(int iB = 0; iB < nBins_e + nBins_mu; iB++){
  	for(int jB = 0; jB < nBins_e + nBins_mu; jB++){
  		chi2 += (FullData[iB]-Prediction[iB])*cov(iB,jB)*(FullData[jB]-Prediction[jB]);
  	}
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    if(!nubar) std::cout << "MBnu Chi2: " << chi2 << std::endl;
    else std::cout << "MBnubar Chi2: " << chi2 << std::endl;
  }
  return chi2;
}
