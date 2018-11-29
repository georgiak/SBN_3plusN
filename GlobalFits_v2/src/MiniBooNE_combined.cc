#include "MiniBooNE_combined.h"
#include "TCanvas.h"


int MiniBooNE_combined::Init(std::string dataLoc, Oscillator osc, bool debug){

  std::string str_data_nue, str_data_numu, str_MC_nue, str_MC_numu, str_fracterrormatrix, str_binboundaries, str_fullosc_nu, str_data_nuebar, str_data_numubar, str_MC_nuebar, str_MC_numubar, str_fullosc_nubar;
  str_binboundaries = dataLoc + "miniboone/miniboone_binboundaries_lowe.txt";
  str_data_nue = dataLoc + "miniboone/data_nue.txt";
  str_data_numu = dataLoc + "miniboone/data_numu.txt";
  str_MC_nue = dataLoc + "miniboone/MC_nue.txt";
  str_MC_numu = dataLoc + "miniboone/MC_numu.txt";
  str_fullosc_nu = dataLoc + "miniboone/numunuefullosc_ntuple.txt";
  str_data_nuebar = dataLoc + "miniboone/miniboone_nuebardata_lowe.txt";
  str_data_numubar = dataLoc + "miniboone/miniboone_numubardata.txt";
  str_MC_nuebar = dataLoc + "miniboone/miniboone_nuebarbgr_lowe.txt";
  str_MC_numubar = dataLoc + "miniboone/miniboone_numubar.txt";
  str_fullosc_nubar = dataLoc + "miniboone/miniboone_nubarfullosc_ntuple.txt";
  str_fracterrormatrix = dataLoc + "miniboone/miniboone_full_fractcovmatrix_combined_lowe.txt";

  float *nu_EnuQE = new float[nBins_e + 1];
	float *nu_FOsc_EnuQE = new float[nFOscEvts_nu];
	float *nu_FOsc_EnuTrue = new float[nFOscEvts_nu];
	float *nu_FOsc_LnuTrue = new float[nFOscEvts_nu];
	float *nu_FOsc_weight = new float[nFOscEvts_nu];
	float *nubar_FOsc_EnuQE = new float[nFOscEvts_nubar];
	float *nubar_FOsc_EnuTrue = new float[nFOscEvts_nubar];
	float *nubar_FOsc_LnuTrue = new float[nFOscEvts_nubar];
	float *nubar_FOsc_weight = new float[nFOscEvts_nubar];

  float *FullData_nu = new float[nBins_e + nBins_mu];
  float *FullData_nubar = new float[nBins_e + nBins_mu];

  ifstream file;
  // Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
  file.open(str_data_nue);
  for(int i = 0; i < nBins_e; i++)
    file >> FullData_nu[i];
  file.close();

  // Get measured Numu ccqe events per enuqe bin
  file.open(str_data_numu);
  for(int i = 0; i < nBins_mu; i++)
    file >> FullData_nu[i + nBins_e];
  file.close();

  // Do the same for antineutrinos
  file.open(str_data_nuebar);
  for(int i = 0; i < nBins_e; i++)
    file >> FullData_nubar[i];
  file.close();
  file.open(str_data_numubar);
  for(int i = 0; i < nBins_mu; i++)
    file >> FullData_nubar[i + nBins_e];
  file.close();

  // Combine fulldata into single  vector
  for(int iB = 0; iB < nBins_mu + nBins_e; iB++){
    FullData[iB] = FullData_nu[iB];
    FullData[iB + nBins_e + nBins_mu] = FullData_nubar[iB];
  }

  // Get predicted nue background events per enuqe bin
  file.open(str_MC_nue);
  for(int i = 0; i < nBins_e; i++)
    file >> Background_nu[i];
  file.close();

  // Get predicted numu ccqe events per enuqe bin
  file.open(str_MC_numu);
  for(int i = 0; i < nBins_mu; i++)
    file >> Background_nu[i + nBins_e];
  file.close();

  // Do the same for antineutrinos
  file.open(str_MC_nuebar);
  for(int i = 0; i < nBins_e; i++)
    file >> Background_nubar[i];
  file.close();
  file.open(str_MC_numubar);
  for(int i = 0; i < nBins_mu; i++)
    file >> Background_nubar[i + nBins_e];
  file.close();

  // Get fractional cov matrix for full numu->nue oscillation events
  file.open(str_fracterrormatrix);
  for(int i = 0; i < 2*nBins_e + 2*nBins_e + 2*nBins_mu; i++)
    for(int j = 0; j < 2*nBins_e + 2*nBins_e + 2*nBins_mu; j++)
      file >> Full_fractCovMatrix[i][j];
  file.close();

  file.open(str_binboundaries);
  for(int i = 0; i < nBins_e + 1; i++)
    file >> nu_EnuQE[i];
  file.close();

  // Get nFullOscEvts for full numu->nue osc events after nue cuts
  file.open(str_fullosc_nu);
  for(int iEvt = 0; iEvt < nFOscEvts_nu; iEvt++){
        file >> nu_FOsc_EnuQE[iEvt];
        file >> nu_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nu_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nu_FOsc_weight[iEvt];    // event weight
  }
  file.close();
  file.open(str_fullosc_nubar);
  for(int iEvt = 0; iEvt < nFOscEvts_nubar; iEvt++){
        file >> nubar_FOsc_EnuQE[iEvt];
        file >> nubar_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nubar_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nubar_FOsc_weight[iEvt];    // event weight
  }
  file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
  float dmmax = 100.;
  float dmmin = 0.01;
  float mstep = TMath::Log10(dmmax/dmmin)/100.f;
  float ETru, LTru, dm2;

  for(int mi = 0; mi < 100; mi++){
    dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
    for(int iB = 0; iB < nBins_e; iB++){
      Lib_sinsq_nu[mi][iB] = 0;
      Lib_sin_nu[mi][iB] = 0;
      Lib_sinsq_nubar[mi][iB] = 0;
      Lib_sin_nubar[mi][iB] = 0;

    }
    for(int iFOsc = 0; iFOsc < nFOscEvts_nu; iFOsc++){   // Loop over full oscillation events
      for(int iB = 0; iB < nBins_e; iB++){    // Loop over energy bins to fill the prediction vector pred

        if(nu_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nu_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){
          ETru = nu_FOsc_EnuTrue[iFOsc];
          LTru = nu_FOsc_LnuTrue[iFOsc];

          Lib_sinsq_nu[mi][iB] += nu_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
          Lib_sin_nu[mi][iB] += nu_FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
        }
      }
    }
    for(int iFOsc = 0; iFOsc < nFOscEvts_nubar; iFOsc++){   // Loop over full oscillation events
      for(int iB = 0; iB < nBins_e; iB++){    // Loop over energy bins to fill the prediction vector pred

        if(nubar_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nubar_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){
          ETru = nubar_FOsc_EnuTrue[iFOsc];
          LTru = nubar_FOsc_LnuTrue[iFOsc];

          Lib_sinsq_nubar[mi][iB] += nubar_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
          Lib_sin_nubar[mi][iB] += nubar_FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
        }
      }
    }
  }

	// Lastly, Initialize best fit signal, with which we will weigh our covariance matrix
  neutrinoModel model_bestFit;
  model_bestFit.zero();
  model_bestFit.Ue[0] = 0.34;  model_bestFit.Um[0] = 0.34;    model_bestFit.mNu[0] = .19;
 	oscContribution oscCont_bestFit_nu = getOscContributionsNueApp(model_bestFit, false, true);
  oscContribution oscCont_bestFit_nubar = getOscContributionsNueApp(model_bestFit, true, true);

	for(int iB = 0; iB < nBins_e; iB++){
		for(int iContribution = 0; iContribution < 6; iContribution++){
    	if(oscCont_bestFit_nubar.dm2[iContribution] == 0) {
        Signal_BestFit_nu[iB] += 0;
        Signal_BestFit_nubar[iB] += 0;
      }
			else{
      	dm2 = floor(TMath::Log10(oscCont_bestFit_nu.dm2[iContribution]/.01)/mstep);
        Signal_BestFit_nu[iB] += oscCont_bestFit_nu.aMuE[iContribution]*Lib_sinsq_nu[dm2][iB] + oscCont_bestFit_nu.aMuE_CPV[iContribution]*Lib_sin_nu[dm2][iB];
        Signal_BestFit_nubar[iB] += oscCont_bestFit_nubar.aMuE[iContribution]*Lib_sinsq_nubar[dm2][iB] + oscCont_bestFit_nubar.aMuE_CPV[iContribution]*Lib_sin_nubar[dm2][iB];
			}
		}
	}

  dof = 2*nBins_e + 2*nBins_mu - 1;

  // Initialize output tree
  chi2Nt = new OutTree("MiniBooNE_combined");

  if(debug){
    std::cout << "MB initialized. Bins: " << dof << std::endl;
  }
  return dof;
}

float MiniBooNE_combined::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  std::array < double, nBins_e + nBins_mu > Prediction_nu, Prediction_nubar;
  std::array < double, 2*nBins_e + 2*nBins_mu > Prediction;

  // Initialize contributions from the oscillation probability
  oscContribution oscCont_nu, oscCont_nubar;
  oscCont_nubar = getOscContributionsNueApp(model, true, true);
  oscCont_nu = getOscContributionsNueApp(model, false, true);

  full_covMatrix.ResizeTo(2*nBins_e + 2*nBins_e + 2*nBins_mu, 2*nBins_e + 2*nBins_e + 2*nBins_mu);
  full_covMatrix.Zero();
  covMatrix.ResizeTo(2*nBins_e + 2*nBins_mu, 2*nBins_e + 2*nBins_mu);
  covMatrix.Zero();

  float mstep = TMath::Log10(100./.01)/float(100);
  int dm2;
  for(int iB = 0; iB < nBins_e; iB++){
    Signal_nu[iB] = 0;
    Signal_nubar[iB] = 0;

  	for(int iContribution = 0; iContribution < 6; iContribution++){
  		if(oscCont_nu.dm2[iContribution] == 0){
        Signal_nu[iB] += 0;
        Signal_nubar[iB] += 0;
      }
  		else{
  			dm2 = floor(TMath::Log10(oscCont_nu.dm2[iContribution]/.01)/mstep);
  			Signal_nu[iB] += oscCont_nu.aMuE[iContribution]*Lib_sinsq_nu[dm2][iB] + oscCont_nu.aMuE_CPV[iContribution]*Lib_sin_nu[dm2][iB];
        Signal_nubar[iB] += oscCont_nubar.aMuE[iContribution]*Lib_sinsq_nubar[dm2][iB] + oscCont_nubar.aMuE_CPV[iContribution]*Lib_sin_nubar[dm2][iB];
  		}
  	}
  }
  // Divide signal prediction by the number of fullosc events
  for(int iB = 0; iB < nBins_e; iB++){
  	Signal_nu[iB] /= double(nFOscEvts_nu);
    Signal_nubar[iB] /= double(nFOscEvts_nubar);
		Signal_BestFit_nu[iB] /= double(nFOscEvts_nu);
    Signal_BestFit_nubar[iB] /= double(nFOscEvts_nubar);
  }

  // Now, scale the fractional cov matrix to our signal and prediction vectors.
  int modesize = nBins_e + nBins_e + nBins_mu;
  int modesize_small = nBins_e + nBins_mu;
  // It's a 6x6 matrix and is therefore an absolute nightmare, but let's have at it!
  for(int iB = 0; iB < 2*nBins_e + 2*nBins_e + 2*nBins_mu; iB++){
    for(int jB = 0; jB < 2*nBins_e + 2*nBins_e + 2*nBins_mu; jB++){
      full_covMatrix(iB,jB) = Full_fractCovMatrix[iB][jB];
      //std::cout << iB << " " << jB << " : " << Full_fractCovMatrix[iB][jB];
      if(iB < nBins_e){
        full_covMatrix(iB,jB) *= Signal_BestFit_nu[iB];
        //std::cout << " * (sig nu) " << Signal_BestFit_nu[iB];
      }
      else if(iB < nBins_e + nBins_e + nBins_mu){
        full_covMatrix(iB,jB) *= Background_nu[iB - nBins_e];
        //std::cout << " * (bkg nu) "  << Background_nu[iB - nBins_e];
      }
      else if(iB < 2*nBins_e + nBins_e + nBins_mu){
        full_covMatrix(iB,jB) *= Signal_BestFit_nubar[iB - modesize];
        //std::cout << " * (sig nubar) " << Signal_BestFit_nubar[iB - modesize];
      }
      else{
        full_covMatrix(iB,jB) *= Background_nubar[iB - modesize - nBins_e];
        //std::cout << " * (bkg nubar)" << Background_nubar[iB - modesize - nBins_e];
      }

      if(jB < nBins_e){
        full_covMatrix(iB,jB) *= Signal_BestFit_nu[jB];
        //std::cout << " * (sig nu) " << Signal_BestFit_nu[jB] << std::endl;
      }
      else if(jB < nBins_e + nBins_e + nBins_mu){
        full_covMatrix(iB,jB) *= Background_nu[jB - nBins_e];
        //std::cout << " * (bkg nu) "  << Background_nu[jB - nBins_e] << std::endl;
      }
      else if(jB < 2*nBins_e + nBins_e + nBins_mu){
        full_covMatrix(iB,jB) *= Signal_BestFit_nubar[jB - modesize];
        //std::cout << " * (sig nubar) " << Signal_BestFit_nubar[jB - modesize] <<  std::endl;
      }
      else{
        full_covMatrix(iB,jB) *= Background_nubar[jB - modesize - nBins_e];
        //std::cout << " * (bkg nubar) " << Background_nubar[jB - modesize - nBins_e] << std::endl;
      }
    }
  }

  // Phew... now go add statistical error
  for(int iB = 0; iB < nBins_e + nBins_e + nBins_mu; iB++){
  	for(int jB = 0; jB < nBins_e + nBins_e + nBins_mu; jB++){
  		if(iB < nBins_e && jB == iB)
        full_covMatrix(iB,iB) += Signal_BestFit_nu[iB];
      if(iB < 2*nBins_e + nBins_e + nBins_mu && iB >= modesize && jB == iB)
        full_covMatrix(iB,iB) += Signal_BestFit_nubar[iB];
    }
  }

  // That was horrible. Now collapse the 6x6 matrix into a 4x4
  for(int iB = 0; iB < nBins_e + nBins_mu; iB++){
  	for(int jB = 0; jB < nBins_e + nBins_mu; jB++){
      if(iB < nBins_e && jB < nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB][jB] + full_covMatrix[iB + nBins_e][jB] + full_covMatrix[iB][jB + nBins_e] + full_covMatrix[iB + nBins_e][jB + nBins_e];
        covMatrix(iB + modesize_small,jB) = full_covMatrix[iB + modesize][jB] + full_covMatrix[iB + modesize + nBins_e][jB] + full_covMatrix[iB + modesize][jB + nBins_e] + full_covMatrix[iB + modesize + nBins_e][jB + nBins_e];
        covMatrix(iB,jB + modesize_small) = full_covMatrix[iB][jB + modesize] + full_covMatrix[iB + nBins_e][jB + modesize] + full_covMatrix[iB][jB + modesize + nBins_e] + full_covMatrix[iB + nBins_e][jB + modesize + nBins_e];
        covMatrix(iB + modesize_small,jB + modesize_small) = full_covMatrix[iB + modesize][jB + modesize] + full_covMatrix[iB + modesize + nBins_e][jB + modesize] + full_covMatrix[iB + modesize][jB + modesize + nBins_e] + full_covMatrix[iB + modesize + nBins_e][jB + modesize + nBins_e];
      }
  		else if(iB < nBins_e && jB >= nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB][jB + nBins_e] + full_covMatrix[iB + nBins_e][jB + nBins_e];
        covMatrix(iB + modesize_small,jB) = full_covMatrix[iB + modesize][jB + nBins_e] + full_covMatrix[iB + modesize + nBins_e][jB + nBins_e];
        covMatrix(iB,jB + modesize_small) = full_covMatrix[iB][jB + modesize + nBins_e] + full_covMatrix[iB + nBins_e][jB + modesize + nBins_e];
        covMatrix(iB + modesize_small,jB + modesize_small) = full_covMatrix[iB + modesize][jB + modesize + nBins_e] + full_covMatrix[iB + modesize + nBins_e][jB + modesize + nBins_e];
  		}
  		else if(iB >= nBins_e && jB < nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB + nBins_e][jB] + full_covMatrix[iB + nBins_e][jB + nBins_e];
        covMatrix(iB + modesize_small,jB) = full_covMatrix[iB + modesize + nBins_e][jB] + full_covMatrix[iB + modesize + nBins_e][jB + nBins_e];
        covMatrix(iB,jB + modesize_small) = full_covMatrix[iB + nBins_e][jB + modesize] + full_covMatrix[iB + nBins_e][jB + modesize + nBins_e];
        covMatrix(iB + modesize_small,jB + modesize_small) = full_covMatrix[iB + modesize + nBins_e][jB + modesize] + full_covMatrix[iB + modesize + nBins_e][jB + modesize + nBins_e];
  		}
  		else if(iB >= nBins_e && jB >= nBins_e){
  			covMatrix(iB,jB) = full_covMatrix[iB + nBins_e][jB + nBins_e];
        covMatrix(iB + modesize_small,jB) = full_covMatrix[iB + modesize + nBins_e][jB + nBins_e];
        covMatrix(iB,jB + modesize_small) = full_covMatrix[iB + nBins_e][jB + modesize + nBins_e];
        covMatrix(iB + modesize_small,jB + modesize_small) = full_covMatrix[iB + modesize + nBins_e][jB + modesize + nBins_e];
  	  }
    }
  }

  // Jesus... now let's invert the covariance matrix
  cov.ResizeTo(2*modesize_small,2*modesize_small);
  cov = covMatrix.Invert();

  for(int iB = 0; iB < modesize_small; iB++){
    Prediction_nu[iB] = Background_nu[iB];
    Prediction_nubar[iB] = Background_nubar[iB];
  }

  for(int iB = 0; iB < nBins_e; iB++){
  	Prediction_nu[iB] += Signal_nu[iB];
    Prediction_nubar[iB] += Signal_nubar[iB];
  }

  // Now combine the nubar and nu into one vector for the prediction and FullData
  for(int iB = 0; iB < nBins_e + nBins_mu; iB++){
    Prediction[iB] = Prediction_nu[iB];
    Prediction[iB + nBins_e + nBins_mu] = Prediction_nubar[iB];
  }

  // Finally, let's put everything together and calculate the chisq
  for(int iB = 0; iB < 2*nBins_e + 2*nBins_mu; iB++){
  	for(int jB = 0; jB < 2*nBins_e + 2*nBins_mu; jB++){
  		chi2 += (FullData[iB]-Prediction[iB])*cov(iB,jB)*(FullData[jB]-Prediction[jB]);
      //std::cout << FullData[iB] << " " << Prediction[iB] << " " << cov(iB,jB) << std::endl;
  	}
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    std::cout << "MiniBooNE Chi2: " << chi2 << std::endl;
  }
  return chi2;
}
