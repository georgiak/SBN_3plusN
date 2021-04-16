#include "MicroBooNE_dis.h"

double mnu_lowbound(.1), mnu_hibound(1000.00);

int MicroBooNE_dis::Init(std::string dataLoc, Oscillator osc, bool debug){

  //////////////////////////////////////////
  shapeonly=true;
  signalInject=true;
  double sin22th_inject(0.0), dm2_inject(1.0);
  dm2_precalc_density = osc.GridSize();
  //////////////////////////////////////////

  Background.resize(nBins);
  FullData.resize(nBins);
  Full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));

  float EnuBinEdges[nBins+1];

  float nu_EnuQE[nMC];
  float nu_EnuTrue[nMC];
  float nu_LnuTrue[nMC];
  float nu_weight[nMC];
  float potweight, fullwgt;

  std::string s_datatag = "apr1_Enu_1m1p";
  std::string s_variable = "Enu_1m1p";
  //std::string s_variable = "Muon_Edep";

  ifstream file;
  file.open(dataLoc+"uboone/"+s_variable+"/microboone_mc_"+s_datatag+".txt");
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << dataLoc+"uboone/microboone_mc_feb16.csv" << std::endl;
	for(short i = 0; i < nMC; i++){
		file >> nu_EnuQE[i];
    file >> nu_EnuTrue[i];
    file >> nu_LnuTrue[i];
    file >> potweight;
    file >> fullwgt;
    nu_weight[i] = potweight*fullwgt;
  }
	file.close();

  file.open(dataLoc+"uboone/"+s_variable+"/microboone_data_histo_"+s_datatag+".txt");
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << dataLoc+"uboone/microboone_data_histo_feb16.csv" << std::endl;
  for(short i = 0; i < nBins; i++){
    file >> FullData[i];
  }
  file.close();

  file.open(dataLoc+"uboone/"+s_variable+"/microboone_bkg_histo_"+s_datatag+".txt");
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << dataLoc+"uboone/microboone_bkg_histo_feb16.csv" << std::endl;
  for(short i = 0; i < nBins; i++){
      file >> Background[i];
  }
  file.close();

  std::cout << "BINEDGES: ";
  file.open(dataLoc+"uboone/"+s_variable+"/microboone_binedges_"+s_datatag+".txt");
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << dataLoc+"uboone/microboone_binedges_feb16.csv" << std::endl;
  for(short i = 0; i <= nBins; i++){
      file >> EnuBinEdges[i];
      std::cout << ", " << EnuBinEdges[i];
  }
  file.close();
  std::cout << std::endl;

	file.open(dataLoc+"uboone/"+s_variable+"/microboone_fracsysmatrix_"+s_datatag+".txt");
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << dataLoc+"uboone/microboone_fracsysmatrix_feb16.csv" << std::endl;
	for(short i = 0; i < nBins; i++)
		for(short j = 0; j < nBins; j++){
			file >> Full_fractCovMatrix[i][j];
    }
	file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
	float mstep = TMath::Log10(pow(mnu_hibound,2)/pow(mnu_lowbound,2))/float(dm2_precalc_density);
	Libdis_sinsq.resize(dm2_precalc_density, std::vector<float>(nBins));
	Libdis_noosc.resize(nBins);

  MCStatSquared_lib.resize(dm2_precalc_density,std::vector<float>(nBins));
  MCStatSquared_noosc.resize(nBins);

	for(int mi = 0; mi < dm2_precalc_density; mi++){
		float dm2 = pow(10,((mi+1.)/float(dm2_precalc_density)*TMath::Log10(pow(mnu_hibound,2)/pow(mnu_lowbound,2)) + TMath::Log10(pow(mnu_lowbound,2))));
    std::cout << "DM2222: " << dm2 << std::endl;
		for(int iB = 0; iB < nBins; iB++){
			Libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				Libdis_noosc[iB] = 0;
		}

		for(int imc = 0; imc < nMC; imc++){   // Loop over mc
      for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	      if(nu_EnuQE[imc] > EnuBinEdges[iB] && nu_EnuQE[imc] < EnuBinEdges[iB+1]){

					float ETru = nu_EnuTrue[imc];
					float LTru = nu_LnuTrue[imc];
					Libdis_sinsq[mi][iB] += nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2); // oscillated events in each bin
          MCStatSquared_lib[mi][iB] += pow(nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2),2);
					if(mi == 0){
						Libdis_noosc[iB] += nu_weight[imc];  // unoscillated events in each bin
            MCStatSquared_noosc[iB] += pow(nu_weight[imc],2);
          }
        }
      }
    }
  }
  dof = nBins;

  // If you want to do a signal injection, I won't stop you. We'll just replace our data.
  if(signalInject){
    std::cout << "INJECTING SIGNAL AT DM2: " << dm2_inject << " SIN22THETA: " << sin22th_inject << std::endl;
    for(int iB = 0; iB < nBins; iB ++){
      std::cout << "DATA: " << FullData[iB] << " ";
      FullData[iB] = Libdis_noosc[iB] + Background[iB];
      std::cout << "NULL: " << FullData[iB] << "  |  " << Libdis_noosc[iB] + Background[iB] << std::endl;
    }
    float mstep = TMath::Log10(pow(mnu_hibound,2)/pow(mnu_lowbound,2))/float(dm2_precalc_density);
    int dm2;
    for(int iB = 0; iB < nBins; iB++){
      dm2 = floor(TMath::Log10(dm2_inject/pow(mnu_lowbound,2))/mstep);
      FullData[iB] -= sin22th_inject*Libdis_sinsq[dm2][iB];
    }
  }

  // Initialize output tree
  chi2Nt = new OutTree("MicroBooNE");

  if(debug){
    std::cout << "MicroBooNE initialized. Bins: " << nBins << std::endl;
  }
  return dof;
}

float MicroBooNE_dis::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  std::vector <  float > Prediction, PredictionNoOsc, MCStatSquared;
  Prediction.resize(nBins);
  PredictionNoOsc.resize(nBins);
  MCStatSquared.resize(nBins);

  // Initialize contributions from osc probability
	double sin22th = model.ProbAmp("mumu");

  covMatrix.ResizeTo(nBins, nBins);
  covMatrix.Zero();

  float minEBins[nBins], maxEBins[nBins];

	for(int iB = 0; iB < nBins; iB ++){
		PredictionNoOsc[iB] = Libdis_noosc[iB] + Background[iB];
	}

	float mstep = TMath::Log10(pow(mnu_hibound,2)/pow(mnu_lowbound,2))/float(dm2_precalc_density);
	int dm2;

	for(int iB = 0; iB < nBins; iB++){
    dm2 = floor(TMath::Log10(model.Dm2()/pow(mnu_lowbound,2))/mstep);
    Prediction[iB] = PredictionNoOsc[iB] - sin22th*Libdis_sinsq[dm2][iB];
    MCStatSquared[iB] = MCStatSquared_noosc[iB] - pow(sin22th,2) * MCStatSquared_lib[dm2][iB];
	}

  if(shapeonly){
    double obsIntegral(0.0), mcIntegral(0.0), covIntegral(0.0), fnorm;
    for(int iB = 0; iB < nBins; iB++){
      obsIntegral += FullData[iB];
      mcIntegral += Prediction[iB];
      for(int jB = 0; jB < nBins;jB++){
        covIntegral += Full_fractCovMatrix[iB][jB];
      }
    }

    fnorm = covIntegral/pow(mcIntegral,2);
    for(int iB = 0; iB < nBins; iB++){
      Prediction[iB] *= (obsIntegral/mcIntegral); // normalize prediction
    }

    for(int iB = 0; iB < nBins; iB++){
  		for(int jB = 0; jB < nBins; jB++){
  			covMatrix[iB][jB] = (Full_fractCovMatrix[iB][jB] - fnorm) * Prediction[iB] * Prediction[jB]; // remove normalization component of cov matrix
  			if(iB == jB){
  				covMatrix[iB][jB] += MCStatSquared[iB];//Prediction[iB];  // Add statistical error of signal prediction
  			}
  		}
  	}
  }
  else{
    for(int iB = 0; iB < nBins; iB++){
		  for(int jB = 0; jB < nBins; jB++){
			  covMatrix[iB][jB] = Full_fractCovMatrix[iB][jB] * Prediction[iB] * Prediction[jB];
			  // Add statistical error of signal prediction
			  if(iB == jB){
				  covMatrix[iB][jB] += MCStatSquared[iB];//Prediction[iB];  // Add statistical error of signal prediction
			  }
		  }
	  }
  }
  //std::cout << std::endl;
	// Now, let's invert this bad boy
	cov.ResizeTo(nBins,nBins);




  std::cout << "err: ";
  for(int iB = 0; iB < nBins; iB++){
		std::cout << sqrt(covMatrix[iB][iB]) << ", ";
	}
  std::cout << std::endl;


  covMatrix[0][0] = .001; // fix null term in cov matrix (this is unused since first bin is empty so it does not matter but it stops me from getting errors)
	cov = covMatrix.Invert();

	// Finally, let's put everything together and calculate the chisq
  double xcheck = 0.0;
	for(int iB = 0; iB < nBins; iB++){
		for(int jB = 0; jB < nBins; jB++){
      xcheck += FullData[iB] - Prediction[iB];
      if(Prediction[iB]>0 && Prediction[jB]>0)
			   chi2 += (FullData[iB] - Prediction[iB]) * cov[iB][jB] * (FullData[jB] - Prediction[jB]);
		}
	}

  std::cout << "XCK: " << xcheck << std::endl;
  /*
  // statsonly
  for(int iB = 0; iB < nBins; iB++){
		if(Prediction[iB]>0)
      chi2 += pow((FullData[iB] - Prediction[iB]),2)/MCStatSquared[iB];
	}
  */


  // Print out spectra:
  std::cout << "data: ";
  for(int iB = 0; iB < nBins; iB++){
		std::cout << FullData[iB] << ", ";
	}
  std::cout << std::endl;

  std::cout << "pred: ";
  for(int iB = 0; iB < nBins; iB++){
		std::cout << Prediction[iB] << ", ";
	}
  std::cout << std::endl;




  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    std::cout << "MicroBooNE Chi2: " << chi2 << std::endl;
  }

  return chi2;
}
