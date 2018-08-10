#include "CDHS.h"

int CDHS::Init(std::string dataLoc, Oscillator osc, bool debug){

  Dm2Vec.resize(601);
  SinSqDeltaGrid_front.resize(dm2VecMaxDim, std::vector<double>(nBins));
  SinSqDeltaGrid_back.resize(dm2VecMaxDim, std::vector<double>(nBins));
  NoOscGrid.resize(nBins, std::vector<double>(2));

  double temp0[] = {0.985,1.006,0.968,1.148,1.000,1.137,1.155,0.887,1.123,0.973,1.039,1.187,1.196,1.006,0.963};
  double temp1[] = {750.,750.,750.,750.,750.,750.,750.,750.,250.,250.,250.,250.,250.,250.,250.};
  double temp2[] = {150.,150.,150.,150.,150.,150.,150.,150.,100.,100.,100.,100.,100.,100.,100.};
  Observed.resize(nBins);
  M_front.resize(nBins);
  M_back.resize(nBins);
  for(int i = 0; i < nBins; i++){
    Observed[i] = temp0[i];
    M_front[i] = temp2[i];
    M_back[i] = temp1[i];
  }

  double sigmaRatio1[] = {0.066,0.055,0.070,0.075,0.087,0.104,0.135,0.200,0.078,0.075,0.085,0.103,0.120,0.132,0.110};
  TMatrix sigmaRatio_inv(nBins,nBins);
  SigmaRatio.ResizeTo(nBins,nBins);
  double sigmaSyst = 0.025;
  for(int i = 0; i < nBins; i++){
    for(int j = 0; j < nBins; j++){
      if(i == j)  sigmaRatio_inv[i][j] = pow(sigmaRatio1[i],2) + pow(sigmaSyst,2);
      else    sigmaRatio_inv[i][j] = pow(sigmaSyst,2);
    }
  }
  SigmaRatio = sigmaRatio_inv.Invert();

  // The cdhs integrals take forever, so i just ran them from the original fortran code and stored them in a txt file, which i read here!
  ifstream file;
  file.open(dataLoc+"cdhs/cdhs_sinsq.txt");
  for(int i = 0; i < nBins; i++){
    for(int k = 0; k < 601; k++){
      file >> Dm2Vec[k];
      file >> SinSqDeltaGrid_front[k][i];
      file >> SinSqDeltaGrid_back[k][i];
      //std::cout << pack.dm2Vec[k] << " " << pack.sinSqDeltaGrid_front[k][i] << " " << pack.sinSqDeltaGrid_back[k][i] << std::endl;
    }
  }

  file.close();
  file.open(dataLoc+"cdhs/cdhs_noosc.txt");
  for(int i = 0; i < nBins; i++){
    file >> NoOscGrid[i][0];
    file >> NoOscGrid[i][1];
  }
  file.close();

  dof = nBins;

  // Iniitialize Output tree
  chi2Nt = new OutTree("CDHS");

  if(debug) std::cout << "CDHS initialized. Bins: " << nBins << std::endl;

  return dof;
}


float CDHS::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  oscContribution oscCont;
  oscCont = getOscContributionsNumuDis(model);

  double sinSqDeltaVec[dm2VecMaxDim];
  double cdhsDm2Vec[601];
  double nFront[nBins], nBack[nBins];
  double nFront_noOsc[nBins], nBack_noOsc[nBins];
  double sinSq, ratio_theor[nBins];

  ROOT::Math::Interpolator dif(dm2VecMaxDim);

  // Front Detector
  for(int iC = 0; iC < nBins; iC++){
    for(int k = 0; k < 601; k++){
      cdhsDm2Vec[k] = Dm2Vec[k];
      sinSqDeltaVec[k] = SinSqDeltaGrid_front[k][iC];
    }

    dif.SetData(601,cdhsDm2Vec,sinSqDeltaVec);

    nFront[iC] = NoOscGrid[iC][0];

    for(int iContribution = 0; iContribution < 6; iContribution++){
      if(oscCont.dm2[iContribution] == 0.)
        sinSq = 0;
      else
        sinSq = dif.Eval(oscCont.dm2[iContribution]);
      nFront[iC] += oscCont.aMuMu[iContribution] * sinSq;
    }
    nFront[iC] *= M_front[iC];
    nFront_noOsc[iC] = M_front[iC] * NoOscGrid[iC][0];
  }

  // Back Detector
  for(int iC = 0; iC < nBins; iC++){
    for(int k = 0; k < 601; k++){
      sinSqDeltaVec[k] = SinSqDeltaGrid_back[k][iC];
    }

    dif.SetData(601,cdhsDm2Vec,sinSqDeltaVec);

    nBack[iC] = NoOscGrid[iC][1];

    for(int iContribution = 0; iContribution < 6; iContribution++){
      if(oscCont.dm2[iContribution] == 0.)
        sinSq = 0;
      else
        sinSq = dif.Eval(oscCont.dm2[iContribution]);
      nBack[iC] += oscCont.aMuMu[iContribution] * sinSq;
    }
    nBack[iC] *= M_back[iC];
    nBack_noOsc[iC] = M_back[iC] * NoOscGrid[iC][1];
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

  if(debug) std::cout << "CDHS Chi2: " << chi2 << std::endl;

  return chi2;
}
