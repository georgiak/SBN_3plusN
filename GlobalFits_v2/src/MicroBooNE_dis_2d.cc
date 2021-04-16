#include "MicroBooNE_dis.h"

int MicroBooNE_dis_2d::Init(std::string dataLoc, Oscillator osc, bool debug){

  //////////////////////////////////////////
  shapeonly=false;
  signalInject=false;
  int dm2_precalc_density = osc.GridSize();
  double sin22th_inject(0.0), dm2_inject(1.0);

  std::string s_datatag = "apr15";
  std::string s_variable1 = "Muon_Edep";
  std::string s_variable2 = "Lepton_CosTheta";
  const int nBins_v1(10);
  const int nBins_v2(3);

  std::cout << "2d Fittin' " << s_variable1 << " " << s_variable2 << std::endl;
  //////////////////////////////////////////
  nBins = nBins_v1*nBins_v2;

  // Instead of full 2d, we're going to flatten to a long, 1d array for easier management. same shit of course

  Background.resize(nBins_v1*nBins_v2);
  FullData.resize(nBins_v1*nBins_v2);

  // stats only  for now

  float binEdges_var1[nBins_v1+1];
  float binEdges_var2[nBins_v2+1];
  std::string s_infile;

  float nu_Var1[nMC];
  float nu_Var2[nMC];
  float nu_EnuTrue[nMC];
  float nu_LnuTrue[nMC];
  float nu_weight[nMC];
  float potweight, fullwgt;

  ifstream file;
  s_infile = dataLoc+"uboone/"+s_variable1+"_"+s_variable2+"_"+s_datatag+"_MC.txt";
  file.open(s_infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_infile << std::endl;
	for(short i = 0; i < nMC; i++){
		file >> nu_Var1[i];
    file >> nu_Var2[i];
    file >> nu_EnuTrue[i];
    file >> nu_LnuTrue[i];
    file >> potweight;
    file >> fullwgt;
    nu_weight[i] = potweight*fullwgt;
  }
	file.close();

  std::cout << "BINEDGES V1: ";
  s_infile = dataLoc+"uboone/"+s_variable1+"_"+s_datatag+"_binedges.txt";
  file.open(s_infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_infile << std::endl;
  for(short i = 0; i <= nBins_v1; i++){
      file >> binEdges_var1[i];
      std::cout << ", " << binEdges_var1[i];
  }
  file.close();
  std::cout << std::endl;

  std::cout << "BINEDGES V2: ";
  s_infile = dataLoc+"uboone/"+s_variable2+"_"+s_datatag+"_binedges.txt";
  file.open(s_infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_infile << std::endl;
  for(short i = 0; i <= nBins_v2; i++){
      file >> binEdges_var2[i];
      std::cout << ", " << binEdges_var2[i];
  }
  file.close();
  std::cout << std::endl;


  s_infile = dataLoc+"uboone/"+s_variable1+"_"+s_variable2+"_"+s_datatag+"_data.txt";
  file.open(s_infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_infile << std::endl;
  for(short i = 0; i < nBins; i++){
    file >> FullData[i];
  }
  file.close();

  s_infile = dataLoc+"uboone/"+s_variable1+"_"+s_variable2+"_"+s_datatag+"_bkg.txt";
  file.open(s_infile);
  if(!file.is_open()) std::cout << "ERROR: BAD FILE NAME: " << s_infile << std::endl;
  for(short i = 0; i < nBins; i++){
    file >> Background[i];
  }
  file.close();

  // To save drastically on the chi2 calculation, precompute all the sines and sine-squareds now!
	float mstep = TMath::Log10(100./.01)/float(dm2_precalc_density);
	Libdis_sinsq.resize(dm2_precalc_density, std::vector<float>(nBins));
	Libdis_noosc.resize(nBins);

  MCStatSquared_lib.resize(dm2_precalc_density,std::vector<float>(nBins));
  MCStatSquared_noosc.resize(nBins);

  int ind, indv1, indv2;
	for(int mi = 0; mi < dm2_precalc_density; mi++){
		float dm2 = pow(10,((mi+.5)/float(dm2_precalc_density)*TMath::Log10(100./.01) + TMath::Log10(.01)));
		for(int iB = 0; iB < nBins; iB++){
			Libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				Libdis_noosc[iB] = 0;
		}

		for(int imc = 0; imc < nMC; imc++){   // Loop over mc
      indv1 = -1;
      indv2 = -1;
      for(int iB = 0; iB < nBins_v1; iB++){    // Loop over bins to fill the prediction vector pred
        if(nu_Var1[imc] > binEdges_var1[iB] && nu_Var1[imc] < binEdges_var1[iB+1]){
          indv1=iB;
          break;
        }
      }
      for(int iB = 0; iB < nBins_v2; iB++){    // Loop over bins to fill the prediction vector pred
        if(nu_Var2[imc] > binEdges_var2[iB] && nu_Var2[imc] < binEdges_var2[iB+1]){
          indv2=iB;
          break;
        }
      }

			float ETru = nu_EnuTrue[imc];
			float LTru = nu_LnuTrue[imc];
      ind = indv1 * nBins_v2 + indv2;

			Libdis_sinsq[mi][ind] += nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2); // oscillated events in each bin
      MCStatSquared_lib[mi][ind] += pow(nu_weight[imc]*pow(sin(1.267*dm2*LTru/ETru),2),2);
			if(mi == 0){
				Libdis_noosc[ind] += nu_weight[imc];  // unoscillated events in each bin
        MCStatSquared_noosc[ind] += pow(nu_weight[imc],2);
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
    float mstep = TMath::Log10(100./.01)/float(100);
    int dm2;
    for(int iB = 0; iB < nBins; iB++){
      dm2 = floor(TMath::Log10(dm2_inject/.01)/mstep);
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

float MicroBooNE_dis_2d::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  std::vector <  float > Prediction, PredictionNoOsc, MCStatSquared;
  Prediction.resize(nBins);
  PredictionNoOsc.resize(nBins);
  MCStatSquared.resize(nBins);

  // Initialize contributions from osc probability
	double sin22th = model.ProbAmp("mumu");
  //model.Print();
  //std::cout << "PROBAMP: " << sin22th << std::endl;

	for(int iB = 0; iB < nBins; iB ++){
		PredictionNoOsc[iB] = Libdis_noosc[iB] + Background[iB];
	}

	float mstep = TMath::Log10(100./.01)/float(100);
	int dm2;

	for(int iB = 0; iB < nBins; iB++){
    dm2 = floor(TMath::Log10(model.Dm2()/.01)/mstep);
    Prediction[iB] = PredictionNoOsc[iB] - sin22th*Libdis_sinsq[dm2][iB];
    MCStatSquared[iB] = MCStatSquared_noosc[iB] - pow(sin22th,2) * MCStatSquared_lib[dm2][iB];
	}

  if(shapeonly){
    double obsIntegral(0.0), mcIntegral(0.0), covIntegral(0.0), fnorm;
    for(int iB = 0; iB < nBins; iB++){
      obsIntegral += FullData[iB];
      mcIntegral += Prediction[iB];
    }
    for(int iB = 0; iB < nBins; iB++){
      Prediction[iB] *= (obsIntegral/mcIntegral); // normalize prediction
    }
  }

  // Finally, let's put everything together and calculate the chisq
	for(int iB = 0; iB < nBins; iB++){
    if(Prediction[iB]>0)
		  chi2 += pow(FullData[iB] - Prediction[iB],2)/MCStatSquared[iB];
	}

  /*
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
  */

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug){
    std::cout << "MicroBooNE Chi2: " << chi2 << std::endl;
  }

  return chi2;
}
