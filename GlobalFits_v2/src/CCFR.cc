#include "CCFR.h"

int CCFR::Init(std::string dataLoc, Oscillator osc, bool debug){

  sinSqDeltaGrid_front.resize(dm2VecMaxDim, std::vector<double>(nBins));
	sinSqDeltaGrid_back.resize(dm2VecMaxDim, std::vector<double>(nBins));
  noOscGrid.resize(nBins, std::vector<double>(2));
  Observed.resize(nBins);
  SigmaRatio.ResizeTo(nBins,nBins);

  for(int i = 0; i < dm2VecMaxDim; i++){
      dm2Vec[i] = pow(10,TMath::Log10(.01) + double(i) / (dm2VecMaxDim-1) * TMath::Log10(100./.01));
  }

	double temp0[] = {0.966,1.014,1.039,0.949,0.988,1.026,1.017,0.963,0.959,1.018,0.990,1.021,1.027,0.971,1.028,1.019,0.938,0.945};
	for(int i = 0; i < nBins; i++)
    Observed[i] = temp0[i];

	TMatrix sigmaRatio_inv(nBins,nBins);
	double NuEnergy[] = {40.76,36.42,96.08,53.62,45.20,128.50,61.01,49.78,151.28,70.51,55.05,176.25,82.86,63.00,209.60,61.06,49.46,150.18};
  double sigmaRatio1[] = {0.034,0.032,0.057,0.025,0.026,0.036,0.032,0.033,0.036,0.025,0.029,0.027,0.030,0.037,0.028,0.032,0.031,0.060};
	double sigmaSyst = 0.015;
	for(int i = 0; i < nBins; i++){
    for(int j = 0; j < nBins; j++){
      if(i == j)  sigmaRatio_inv[i][j] = pow(sigmaRatio1[i],2) + pow(sigmaSyst,2);
      else    sigmaRatio_inv[i][j] = pow(sigmaSyst,2);
		}
	}
	SigmaRatio = sigmaRatio_inv.Invert();

	double lBack = 1.116;                      double lFront = 0.715;
	double deltaLBack = .357;                  double deltaLFront = .3565;
	double l1Back = lBack - deltaLBack/2;      double l2Back = lBack + deltaLBack/2;
	double l1Front = lFront - deltaLFront/2;   double l2Front = lFront + deltaLFront/2;

	// Front Detector
	for(int iC = 0; iC < nBins; iC++){
    noOscGrid[iC][0] = (l2Front - l1Front) / (l1Front * l2Front);
	}
	for(int k = 0; k < dm2VecMaxDim; k ++){
		for(int iC = 0; iC < nBins; iC++){
      double Enu = NuEnergy[iC];
		  double delta = 1.27 * dm2Vec[k] / Enu;
		  double num1 = 1. - cos(2. * delta * l1Front) - 2. * delta * l1Front * sineInt(2. * delta * l1Front);
		  double num2 = 1. - cos(2. * delta * l2Front) - 2. * delta * l2Front * sineInt(2. * delta * l2Front);
		  double den1 = 2 * l1Front;  double den2 = 2* l2Front;
		  double term1 = num1 / den1; double term2 = num2 / den2;
		  sinSqDeltaGrid_front[k][iC] = term1 - term2;
    }
	}
	// Back Detector
	for(int iC = 0; iC < nBins; iC++){
	  noOscGrid[iC][1] = (l2Back - l1Back) / (l1Back * l2Back);
	}
	for(int k = 0; k < dm2VecMaxDim; k ++){
    for(int iC = 0; iC < nBins; iC++){
      double Enu = NuEnergy[iC];
	    double delta = 1.27 * dm2Vec[k] / Enu;
	    double num1 = 1. - cos(2. * delta * l1Back) - 2. * delta * l1Back * sineInt(2. * delta * l1Back);
			double num2 = 1. - cos(2. * delta * l2Back) - 2. * delta * l2Back * sineInt(2. * delta * l2Back);
	    double den1 = 2 * l1Back;  double den2 = 2* l2Back;
	    double term1 = num1 / den1; double term2 = num2 / den2;
	    sinSqDeltaGrid_back[k][iC] = term1 - term2;
    }
  }

  dof = nBins;

  //Initialize output tree
  chi2Nt = new OutTree("CCFR");

  if(debug) std::cout << "CCFR initialized. Bins: " << nBins << std::endl;

  return dof;
}

float CCFR::Chi2(Oscillator osc, neutrinoModel model, bool debug){

  float chi2 = 0.f;
  std::vector < float > Prediction;
  Prediction.resize(nBins);

  // Initialize contributions from the osc probability
  oscContribution oscCont;
  oscCont = getOscContributionsNumuDis(model);

  ROOT::Math::Interpolator dif(dm2VecMaxDim);
  double sinSqDeltaVec[dm2VecMaxDim];
  double nFront[nBins], nBack[nBins];
  double nFront_noOsc[nBins], nBack_noOsc[nBins];
  double sinSq, ratio_theor[nBins];

  // Front Detector
  for(int iC = 0; iC < nBins; iC++){
    for(int k = 0; k < dm2VecMaxDim; k++){
      sinSqDeltaVec[k] = sinSqDeltaGrid_front[k][iC];
    }
    dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);

    nFront[iC] = noOscGrid[iC][0];

    for(int iContribution = 0; iContribution < 6; iContribution++){
      if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
      else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
      nFront[iC] += oscCont.aMuMu[iContribution] * sinSq;
    }
    nFront[iC] *= m_front;
    nFront_noOsc[iC] = m_front * noOscGrid[iC][0];
  }

  // Back Detector
  for(int iC = 0; iC < nBins; iC++){
    for(int k = 0; k < dm2VecMaxDim; k++){
      sinSqDeltaVec[k] = sinSqDeltaGrid_back[k][iC];
    }
    dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);

    nBack[iC] = noOscGrid[iC][1];

    for(int iContribution = 0; iContribution < 6; iContribution++){
      if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
      else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
      nBack[iC] += oscCont.aMuMu[iContribution] * sinSq;
    }
    nBack[iC] *= m_back;
    nBack_noOsc[iC] = m_back * noOscGrid[iC][1];
  }

  // Get the predicted ratio
  for(int iC = 0; iC < nBins; iC++){
    ratio_theor[iC] = (nBack[iC]/nFront[iC])/(nBack_noOsc[iC]/nFront_noOsc[iC]);
  }

  // Calculate the chisq
  for(int iC = 0; iC < nBins; iC++){
    for(int jC = 0; jC < nBins; jC++){
      chi2 += (Observed[iC] - ratio_theor[iC])*SigmaRatio[iC][jC]*(Observed[jC] - ratio_theor[jC]);
    }
  }

  // Fill output tree
  chi2Nt->Fill(chi2,dof,model);

  if(debug)
    std::cout << "CCFR Chi2: " << chi2 << std::endl;
  return chi2;
}
