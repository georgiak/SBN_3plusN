#include "MINOS.h"

int MINOS::Init(std::string dataLoc, Oscillator osc, bool debug){

  EnuQE.resize(nBins+1);
  NumubarData.resize(nBins);
  NumubarBkg.resize(nBins);
  fracError.resize(nBins);
  dataErr.resize(nBins);
  EnuQE_ws.resize(nBins_ws+1);
  NumubarData_ws.resize(nBins_ws);
  NumubarBkg_ws.resize(nBins_ws);
  fracError_ws.resize(nBins_ws);
  dataErr_ws.resize(nBins_ws);

  Signal.resize(nBins);
  Prediction.resize(nBins);
  TotalError.resize(nBins);

  // Minos looks easy, but DON'T BE FOOLED.
  // Right sign
  double temp0[] = {0, 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 15000., 20000.};
  double temp1[] = {0.14, 6.0, 13.1, 41.0, 29.0, 19.1, 7.9, 11.1, 11.1, 4.0, 4.8, 4.4};
  double temp2[] = {1.2, 20.9, 50.5, 56.8, 34.5, 17.8, 11.8, 9.1, 7.8, 7.0, 5.0, 2.7};
  double temp3[] = {0.020, 0.020, 0.06, 0.05, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.020, 0.020};
  double temp4[] = {1.04, 3.21, 4.30, 7.02, 6.10, 5.12, 3.54, 4.03, 4.03, 2.78, 1.31, 1.25};
  for(int i = 0; i < nBins; i++){
    EnuQE[i] = temp0[i];
    NumubarData[i] = temp1[i];
    NumubarBkg[i] = temp2[i];
    fracError[i] = temp3[i];
    dataErr[i] = temp4[i];
  }
  EnuQE[nBins] = temp0[nBins];

  // Wrong sign
  double temp5[] = {0, 2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000., 25000., 30000., 50000.};
  double temp6[] = {0.0, 6.99, 22.00, 16.01, 14.95, 14.94, 13.94, 9.99, 7.94, 3.89, 1.94, 1.93, 0.73};
  double temp7[] = {3.4, 12.72, 16.31, 18.15, 18.99, 16.44, 13.94, 11.14, 8.84, 6.93, 4.38, 2.08, 0.53};
  double temp8[] = {.17, .13, .11, .1, .09, .09, .09, .11, .1, .11, .09, .16, .3};
  double temp9[] = {1.93, 3.29, 5.36, 4.67, 4.54, 4.52, 4.42, 3.79, 3.52, 2.67, 1.25, 1.2, .52};
  for(int i = 0; i < nBins_ws; i++){
    EnuQE_ws[i] = temp5[i];
    NumubarData_ws[i] = temp6[i];
    NumubarBkg_ws[i] = temp7[i];
    fracError_ws[i] = temp8[i];
    dataErr_ws[i] = temp9[i];
  }
  EnuQE_ws[nBins_ws] = temp5[nBins_ws];

  dof = nBins;

  // Initalize our output tree
  chi2Nt = new OutTree("MINOS");

	if(debug) std::cout << "MINOS initialized. Bins: " << nBins << std::endl;

  return dof;
}

float MINOS::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  double minEBins[13]; double maxEBins[13];
  double ETrue[13];
  double LTrue = 735000.;     double LTrue_near = 1040.;

  model.difference();
  double sin41Sq,sin51Sq,sin31Sq,sin53Sq,sin43Sq,sin54Sq,sin61Sq,sin63Sq,sin64Sq,sin65Sq;
  double dm2Atm = 0.00232;    double dm31Sq=dm2Atm;
  double Um3 = 0.71;
  double sin22Theta_atm_err = 0.1;
  double dm43Sq = model.dm41Sq - dm31Sq;
  double dm53Sq = model.dm51Sq - dm31Sq;
  double dm63Sq = model.dm61Sq - dm31Sq;

  // Minos Nubar Right Sign =========
  Signal.assign(nBins,0.);
  Prediction.assign(nBins,0.);
  TotalError.assign(nBins,0.);

  for(int iM = 0; iM < nBins; iM++){
    minEBins[iM] = EnuQE[iM];
    maxEBins[iM] = EnuQE[iM + 1];

    ETrue[iM] = (maxEBins[iM] - minEBins[iM])/2. + minEBins[iM];
    double binWidth = (maxEBins[iM] - minEBins[iM])/1000.;
    double EnuGeV = ETrue[iM]/1000;

    // Check that the oscillations aren't too fast for our resolution
    double nMax = 1.27 * model.dm41Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    double EnuNext = 1.27 * model.dm41Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth){
      sin41Sq = pow(sin(1.27 * model.dm41Sq * LTrue / ETrue[iM]),2);
      sin43Sq = pow(sin(1.27 * dm43Sq * LTrue / ETrue[iM]),2);
    }
    else{
      sin41Sq = 0.5;
      sin43Sq = 0.5;
    }

    nMax = 1.27 * model.dm51Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    EnuNext = 1.27 * model.dm51Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth){
      sin51Sq = pow(sin(1.27 * model.dm51Sq * LTrue / ETrue[iM]),2);
      sin53Sq = pow(sin(1.27 * dm53Sq * LTrue / ETrue[iM]),2);
    }
    else{
      sin51Sq = 0.5;
      sin53Sq = 0.5;
    }

    nMax = 1.27 * model.dm54Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    EnuNext = 1.27 * model.dm54Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin54Sq = pow(sin(1.27 * model.dm54Sq * LTrue / ETrue[iM]),2);
    else
      sin54Sq = 0.5;

    nMax = 1.27 * model.dm65Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    EnuNext = 1.27 * model.dm65Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin65Sq = pow(sin(1.27 * model.dm65Sq * LTrue / ETrue[iM]),2);
    else
      sin65Sq = 0.5;

    nMax = 1.27 * model.dm64Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    EnuNext = 1.27 * model.dm64Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin64Sq = pow(sin(1.27 * model.dm64Sq * LTrue / ETrue[iM]),2);
    else
      sin64Sq = 0.5;

    nMax = 1.27 * model.dm61Sq * (LTrue / 1000.) / (EnuGeV * TMath::Pi());
    EnuNext = 1.27 * model.dm61Sq * (LTrue/1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth){
      sin61Sq = pow(sin(1.27 * model.dm61Sq * LTrue / ETrue[iM]),2);
      sin63Sq = pow(sin(1.27 * dm63Sq * LTrue / ETrue[iM]),2);
    }
    else{
      sin61Sq = 0.5;
      sin63Sq = 0.5;
    }

    sin31Sq = pow(sin(1.27 * dm31Sq * LTrue / ETrue[iM]),2);

    Prediction[iM] = NumubarBkg[iM] * (1 - 4 * ((1 - pow(Um3,2) - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(Um3,2) * sin31Sq + pow(model.Um[0],2) * sin41Sq + pow(model.Um[1],2) * sin51Sq + pow(model.Um[2],2) * sin61Sq)
              + pow(model.Um[0],2) * pow(Um3,2) * sin43Sq + pow(model.Um[1],2) * pow(model.Um[0],2) * sin54Sq + pow(model.Um[1],2) * pow(Um3,2) * sin53Sq
              + pow(model.Um[2],2) * pow(model.Um[1],2) * sin65Sq + pow(model.Um[2],2) * pow(model.Um[0],2) * sin64Sq + pow(model.Um[2],2) * pow(Um3,2) * sin63Sq));

    // Now, do the same horrible thing for the near detector;
    nMax = 1.27 * model.dm41Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm41Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin41Sq = pow(sin(1.27 * model.dm41Sq * LTrue_near / ETrue[iM]),2);
    else
      sin41Sq = 0.5;

    nMax = 1.27 * model.dm51Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm51Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin51Sq = pow(sin(1.27 * model.dm51Sq * LTrue_near / ETrue[iM]),2);
    else
      sin51Sq = 0.5;

    nMax = 1.27 * model.dm54Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm54Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin54Sq = pow(sin(1.27 * model.dm54Sq * LTrue_near / ETrue[iM]),2);
    else
      sin54Sq = 0.5;

    nMax = 1.27 * model.dm65Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm65Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin65Sq = pow(sin(1.27 * model.dm65Sq * LTrue_near / ETrue[iM]),2);
    else
      sin65Sq = 0.5;

    nMax = 1.27 * model.dm64Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm64Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin64Sq = pow(sin(1.27 * model.dm64Sq * LTrue_near / ETrue[iM]),2);
    else
      sin64Sq = 0.5;

    nMax = 1.27 * model.dm61Sq * (LTrue_near / 1000.) / (EnuGeV*TMath::Pi());
    EnuNext = 1.27 * model.dm61Sq * (LTrue_near / 1000.) / ((nMax + 0.5) * TMath::Pi());
    if(abs(EnuGeV - EnuNext) >= binWidth)
      sin61Sq = pow(sin(1.27 * model.dm61Sq * LTrue_near / ETrue[iM]),2);
    else
      sin61Sq = 0.5;

    Signal[iM] = (1 - 4* ((1 - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) * sin41Sq + pow(model.Um[1],2) * sin51Sq + pow(model.Um[2],2) * sin61Sq)
              + pow(model.Um[1],2) * pow(model.Um[0],2) * sin54Sq + pow(model.Um[2],2) * pow(model.Um[1],2) * sin65Sq + pow(model.Um[2],2) * pow(model.Um[0],2) * sin64Sq));

    Prediction[iM] /= Signal[iM];

    if(Prediction[iM] < 0)  Prediction[iM] = 0;

    // Now, fill up the total error while we're still loopin' around
    double error_from_atm_mixing = NumubarBkg[iM] * sin22Theta_atm_err * pow(sin(1.27 * dm2Atm * LTrue / ETrue[iM]),2);
    TotalError[iM] = pow(dataErr[iM],2) + pow(fracError[iM] * Prediction[iM],2) + pow(error_from_atm_mixing,2);
  }

  // Now, calculate the right-sign chi2
  for(int iM = 0; iM < nBins; iM++){
    chi2 += (NumubarData[iM] - Prediction[iM]) * (NumubarData[iM] - Prediction[iM]) / TotalError[iM];
  }

  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "MINOS Chi2: " << chi2 << std::endl;

  return chi2;
}
