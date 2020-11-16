#include "NEOS.h"

/*
Aug 29, Davio Cianci

Fit method: https://arxiv.org/pdf/1709.04294.pdf
Energy Resolution: https://arxiv.org/pdf/1609.03910.pdf
NEOS results: https://arxiv.org/pdf/1610.05134.pdf
Antineutrino fluxes from https://arxiv.org/pdf/1101.2663.pdf and https://arxiv.org/pdf/1106.0687.pdf
IBD XSec from https://arxiv.org/pdf/hep-ph/9903554.pdf

*/


int NEOS::Init(std::string dataLoc, Oscillator osc, bool debug){

  double en, len, prob, norm;
  Cov.ResizeTo(60,60);


  // neos/dayabay data
  Observed = {0.97774, 1.03994, 1.03064, 0.96082, 0.95823, 0.96296, 0.96341, 0.97332, 0.98857,
              0.96860, 0.98643, 0.99329, 1.00274, 0.98689, 1.01509, 1.01616, 0.98262, 0.99360,
              0.99360, 1.00290, 1.01494, 1.03430, 1.02287, 1.00198, 1.00259, 1.00000, 0.99512,
              1.01890, 1.01280, 0.98598, 0.99405, 0.99619, 1.00442, 1.01402, 1.01479, 1.00473,
              1.03171, 1.02576, 1.01128, 1.04253, 1.01280, 1.02454, 1.00869, 1.00488, 0.97561,
              0.99070, 0.97820, 0.96540, 1.03887, 1.02256, 0.97866, 1.01387, 1.00244, 1.05046,
              1.03979, 1.03155, 1.01692, 1.03750, 0.99741, 1.05168};

  // Fluxes
  // NOTE to account for high and low energy, I incorrectly assume linearity before first point and asymptote to zero at the end in order to
  // extend the function across my desired energy range.
  // This will only affect neutrino energies < 2MeV and > 8MeV, which only come in through smearing anyways.
  // I interpolated these in python to make things easier
  ifstream file;
  file.open(dataLoc+"neos_fluxarray.txt");
  for(int iIso = 0; iIso < 4; iIso++){
    for(int iEn = 0; iEn < 200; iEn++){
      file >> Flux_iso[iIso][iEn];
    }
  }
  file.close();

  // For xsec, energy and 10^42 * the IBD cross section
  double xsec_en[] = {0,1.6453, 1.8179, 1.9655, 2.1132, 2.2692, 2.5513, 2.8504, 3.1239, 3.2863, 3.4421, 3.5898, 3.7374, 3.8851, 4.0327, 4.1804, 4.3280, 4.4757,
                      4.6234, 4.7710, 4.9187, 5.0663, 5.2140, 5.3616, 5.5093, 5.6570, 5.8046, 5.9523, 6.0999, 6.2476, 6.3952, 6.5429, 6.8382, 6.9859, 7.1335,
                      7.2812, 7.4288, 7.7242, 7.8718, 8.0195, 8.3718, 8.7578, 8.9054, 9.0531, 9.2007, 9.3484, 9.4961, 9.6213, 9.7914, 9.9122, 12.094};
  double xsec[] = {   0,0.0083717, 0.0064888, 0.021635, 0.051349, 0.083366, 0.13609, 0.21968, 0.31339, 0.37576, 0.43451, 0.48945, 0.56342, 0.63008, 0.70259, 0.77949,
                      0.86079, 0.95234, 1.0410, 1.1223, 1.2304, 1.3259, 1.4189, 1.5470, 1.6474, 1.7565, 1.8920, 2.0050, 2.1258, 2.2729, 2.3943, 2.5259, 2.8227, 2.9508,
                      3.1142, 3.2365, 3.4083, 3.7884, 3.9329, 4.1248, 4.4941, 4.9959, 5.1632, 5.3587, 5.5571, 5.7600, 5.9672, 6.1411, 6.3831, 6.5668, 10.514};
  ROOT::Math::Interpolator dif1(51);
  dif1.SetData(51,xsec_en,xsec);
  for(int i = 0; i < 200; i++){
    en = .05 * i + .05/2.;
    XSec[i] = dif1.Eval(en+.8) * 1e-42;
  }

  // Isotope fission fraction (array in order u235 u238 pu239 pu241) for each detector
  DB_f_iso = {{ { .564, .564, .557, .552},
                { .076, .076, .076, .076},
                { .303, .303, .312, .315},
                { .056, .056, .055, .057}}};
  NEOS_f_iso = { 0.655, 0.072, 0.235, 0.038};

  // Efficiency
  double eff_mu[] = { .8255, .8221, .8573, .8571};
  double eff_mult[] = {.9744, .9747, .9757, .9757};
  DB_eff = {eff_mu[0] * eff_mult[0], eff_mu[1] * eff_mult[1], eff_mu[2] * eff_mult[2], eff_mu[3] * eff_mult[3]};

  // Length  (array for each detector in order AD1 AD2 AD3 AD8) for each reactor at Daya Bay
  DB_l_d = {{ {362.38, 371.76, 903.47, 817.16, 1353.62, 1265.32},
              {357.94, 368.41, 903.35, 816.90, 1354.23, 1265.89},
              {1332.48, 1358.15, 467.57, 489.58, 557.58, 499.2},
              {1337.43, 1362.88, 472.97, 495.35, 558.71, 501.07}}};

  // Active detector mass for Daya Bay
  DB_mass = { 907.f, 916.f, 915.f,  950.f };

  // Smearing matrix from Daya Bay paper
  file.open(dataLoc+"DayaBay_DetectorModel_Etrue_to_Erec.txt");
  for(short i = 0; i < 240; i++){
    for(short j = 0; j < 240; j++){
  	  file >> DB_smearingmatrix[i][j];
    }
  }
  file.close();

  // Generate the smearing matrix for NEOS_H
  // NOTE: this does not account for positrons escaping the detector, which would result in a flat, non-zero baseline
  // to the left of the gaussian peak ( sorta like -^\_ ). Described in one of the papers above. I'll look into it.
  int nobs = 100000;
  double resol;
  std::string st_gaus;

  for(int i = 0; i < 200; i++){
    en = .05 * i + .05/2.;
    resol = 0.05 * sqrt(en) + .12;
    st_gaus = "TMath::Gaus(x," + to_string(en) +","+ to_string(resol) + ")";
    TF1* mygaus = new TF1("mygaus",st_gaus.c_str(),0,10);
    TH1D* shist = new TH1D("shist","gaus",200,0,10);
    shist->FillRandom("mygaus",nobs);
    norm = shist->Integral();
    for(int j = 0; j < 200; j++){
      NEOS_smearingmatrix[j][i] = shist->GetBinContent(j+1)/norm;
    }
    delete shist;
    delete mygaus;
  }

  // Now actually calculate event rates!
  //
  // Daya Bay 3v
  //
  std::array < std::array < double, 200 >, 4 > DB_Predicted_3n_Unsmeared, DB_Predicted_3n_Smeared;
  for(int iDet = 0; iDet < 4; iDet ++){ // Loop through detectors
    for(int iEn = 0; iEn < 200; iEn++){ // Loop through positron energy from 0-10MeV
      DB_Predicted_3n_Unsmeared[iDet][iEn] = 0;

      // Average across each bin to catch fast oscillations (particularly in low bins)
      for(int i = 0; i < binAvg; i++){
        en = (.05 * iEn) + .05 * (i+1)/float(lenAvg+1);

        for(int iReac = 0; iReac < 6; iReac ++){ // Loop through reactors to get osc length
          for(int iIso = 0; iIso < 4; iIso ++){ // Loop through isotopes
            prob = 1.;
            DB_Predicted_3n_Unsmeared[iDet][iEn] += Flux_iso[iIso][iEn] * XSec[iEn] * DB_f_iso[iIso][iDet] * prob * DB_eff[iDet] / pow(DB_l_d[iDet][iReac],2);
          }
        }
      }
      DB_Predicted_3n_Unsmeared[iDet][iEn] /= float(binAvg);
    }
  }
  // Apply energy smearing
  for(int i = 0; i < 200; i++){
    for(int iDet = 0; iDet < 4; iDet++){
      DB_Predicted_3n_Smeared[iDet][i] = 0;
      for(int j = 0; j < 200; j++){
        if(iDet < 2)  norm = DB_norm_EH1;
        else norm = DB_norm_EH2;
        DB_Predicted_3n_Smeared[iDet][i] += DB_Predicted_3n_Unsmeared[iDet][j] * DB_smearingmatrix[i][j] * norm;
      }
    }
  }
  // Lastly, average to appropriate bins and weigh detector contributions by their masses
  for(int iB = 0; iB < nBins; iB++){
    DB_Predicted_3n[iB] = 0;
    for(int iDet = 0; iDet < 4; iDet++){
      // Need to average between two because of the binning of the smearing matrix
      DB_Predicted_3n[iB] += ((DB_Predicted_3n_Smeared[iDet][20+iB*2] + DB_Predicted_3n_Smeared[iDet][20+iB*2 + 1])/2.0) * (DB_mass[iDet] / 922.f);
    }
  }

  //
  // NEOS 3v hypothesis
  //
  std::array < double, 200 > NEOS_Predicted_3n_Unsmeared, NEOS_Predicted_3n_Smeared;
  for(int iEn = 0; iEn < 200; iEn++){           // Loop through energy from 0-10MeV
    NEOS_Predicted_3n_Unsmeared[iEn] = 0;

    // Average across each bin to catch fast oscillations (particularly in low bins)
    for(int i = 0; i < binAvg; i++){
      en = (.05 * iEn) + .05 * (i+1)/float(lenAvg+1);

      // Average across length to account for physical detector size
      for(int j = 0; j < lenAvg; j++){
        len = osc.RanGen.Gaus(24.0,1.5/3.);

        for(int iIso = 0; iIso < 4; iIso ++){ // Loop through isotopes
          prob = 1.;
          NEOS_Predicted_3n_Unsmeared[iEn] += Flux_iso[iIso][iEn] * XSec[iEn] * NEOS_f_iso[iIso] * prob / pow(len,2);
        }
      }
    }
  }
  // Apply energy smearing
  for(int i = 0; i < 200; i++){
    NEOS_Predicted_3n_Smeared[i] = 0;
    for(int j = 0; j < 200; j++){
      NEOS_Predicted_3n_Smeared[i] += NEOS_Predicted_3n_Unsmeared[j] * NEOS_smearingmatrix[i][j];
    }
  }
  // Lastly, average to appropriate bins and get NEOS ratio
  for(int iB = 0; iB < nBins; iB++){
    NEOS_Predicted_3n[iB] = (NEOS_Predicted_3n_Smeared[20+iB*2] + NEOS_Predicted_3n_Smeared[20+iB*2 + 1])/2.0;
  }

  // load up fractional covariance Matrix
  file.open(dataLoc+"neos_cov.txt");
  for(short i = 0; i < 60; i++){
    for(short j = 0; j < 60; j++){
  	  file >> Cov[i][j];
    }
  }
  file.close();
/*
  // Make full cov matrix
  for(int iB = 0; iB < 60; iB++){
    for(int jB = 0; jB < 60; jB++){
      Cov[iB][jB] *= sqrt(StatsError[iB]) * sqrt(StatsError[jB]);
    }
  }

  // Add stats error to cov Matrix
  for(int iB = 0; iB < 60; iB++){
    Cov[iB][iB] += StatsError[iB];
  }
  */

  Cov = Cov.Invert();

  dof = nBins;

  // Iniitalize our output tree
  chi2Nt = new OutTree("NEOS");

	if(debug) std::cout << "NEOS initialized. Bins: " << nBins << std::endl;

  return dof;
}

float NEOS::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  double sin22th = model.ProbAmp("ee");
  double dm2 = model.Dm2();

  std::array < double, nBins > DB_Predicted_4n, DB_Ratio, NEOS_Predicted_4n, NEOS_Ratio;
  double prob, en, len, norm;

  //
  // Daya Bay 4v
  //
  std::array < std::array < double, 200 >, 4 > DB_Predicted_4n_Unsmeared, DB_Predicted_4n_Smeared;
  for(int iDet = 0; iDet < 4; iDet ++){ // Loop through detectors
    for(int iEn = 0; iEn < 200; iEn++){
      DB_Predicted_4n_Unsmeared[iDet][iEn] = 0;

      // Average across each bin to catch fast oscillations (particularly in low bins)
      for(int i = 0; i < binAvg; i++){
        en = (.05 * iEn) + .05 * (i+1)/float(lenAvg+1);

        for(int iReac = 0; iReac < 6; iReac ++){ // Loop through reactors
          for(int iIso = 0; iIso < 4; iIso ++){ // Loop through isotopes
            prob = 1. - sin22th * pow(sin(1.267 * dm2 * DB_l_d[iDet][iReac] / (en+.8)),2);
            DB_Predicted_4n_Unsmeared[iDet][iEn] += Flux_iso[iIso][iEn] * XSec[iEn] * DB_f_iso[iIso][iDet] * prob * DB_eff[iDet] / pow(DB_l_d[iDet][iReac],2);
          }
        }
      }
      DB_Predicted_4n_Unsmeared[iDet][iEn] /= float(binAvg);
    }
  }
  // Apply energy smearing
  for(int i = 0; i < 200; i++){
    for(int iDet = 0; iDet < 4; iDet++){
      DB_Predicted_4n_Smeared[iDet][i] = 0;
      for(int j = 0; j < 200; j++){
        if(iDet < 2) norm = DB_norm_EH1;
        else  norm = DB_norm_EH2;
        DB_Predicted_4n_Smeared[iDet][i] += DB_Predicted_4n_Unsmeared[iDet][j] * DB_smearingmatrix[i][j] * norm;
      }
    }
  }
  // Lastly, average to appropriate bins and get DayaBay ratio
  for(int iB = 0; iB < nBins; iB++){
    DB_Predicted_4n[iB] = 0;
    for(int iDet = 0; iDet < 4; iDet++){
      DB_Predicted_4n[iB] += ((DB_Predicted_4n_Smeared[iDet][20+iB*2] + DB_Predicted_4n_Smeared[iDet][20+iB*2 + 1])/2.0) * (DB_mass[iDet] / 922.f);
    }
    // Get ratio
    DB_Ratio[iB] = DB_Predicted_4n[iB]/DB_Predicted_3n[iB];
    //std::cout << "DB final: " << DB_Predicted_4n[iB] << " "  <<  NEOS_Predicted_3n[iB] << std::endl;
  }


  //
  // NEOS 4v
  //
  std::array < double, 200 > NEOS_Predicted_4n_Unsmeared, NEOS_Predicted_4n_Smeared;
  for(int iEn = 0; iEn < 200; iEn++){
    NEOS_Predicted_4n_Unsmeared[iEn] = 0;

    // Avg across each bin to catch fast oscillations
    for(int i = 0; i < binAvg; i++){
      en = (.05 * iEn) + .05 * (i+1)/float(lenAvg+1);

      // Avg across length to account for detector
      for(int j = 0; j < lenAvg; j++){
        len = osc.RanGen.Gaus(24.0,1.5/3.);

        for(int iIso = 0; iIso < 4; iIso ++){ // Loop through isotopes
          prob = 1. - sin22th * pow(sin(1.267 * dm2 * len / (en+.8)),2);
          NEOS_Predicted_4n_Unsmeared[iEn] += Flux_iso[iIso][iEn] * XSec[iEn] * NEOS_f_iso[iIso] * prob / pow(len,2);
        }
      }
    }
  }
  // Apply energy smearing
  for(int i = 0; i < 200; i++){
    NEOS_Predicted_4n_Smeared[i] = 0;
    for(int j = 0; j < 200; j++){
      NEOS_Predicted_4n_Smeared[i] += NEOS_Predicted_4n_Unsmeared[j] * NEOS_smearingmatrix[i][j];
    }
  }
  // Lastly, average to appropriate bins and get NEOS ratio
  for(int iB = 0; iB < nBins; iB++){
    NEOS_Predicted_4n[iB] = (NEOS_Predicted_4n_Smeared[20+iB*2] + NEOS_Predicted_4n_Smeared[20+iB*2 + 1])/2.0;

    // Get ratio
    NEOS_Ratio[iB] = NEOS_Predicted_4n[iB]/NEOS_Predicted_3n[iB];
    //std::cout << "NEOS final: " << NEOS_Predicted_4n[iB] << " "  <<  NEOS_Predicted_3n[iB] << std::endl;
  }

  //
  // NOW WE HAVE THE FINAL RATIO FINALLY! let's see how off we are.
  //
  std::array < double, nBins > Ratio;
  for(int iB = 0; iB < nBins; iB++){
    Ratio[iB] = NEOS_Ratio[iB]/DB_Ratio[iB];
    //std::cout << "RATIO: " << Ratio[iB] << std::endl;
  }

  // Calc chi2
  for(int iB = 0; iB < nBins; iB++){
    for(int jB = 0; jB < nBins; jB++){
      chi2 += (Observed[iB] - Ratio[iB]) * Cov[iB][jB] * (Observed[jB] - Ratio[jB]);
      //if(chi2 < 0)
        //std:cout << "NEG: " << (Observed[iB] - Ratio[iB]) << " " <<  (Observed[jB] - Ratio[jB]) << " COV: " << Cov[iB][jB] << std::endl;
    }
  }

/*
  std::array<double,60> Cov_diag = {
            9.20600000e-06, 9.20600000e-06, 9.20600000e-06, 1.09069839e-05,2.31834075e-05, 2.57525221e-05,
            2.01124903e-05, 1.57774130e-05,1.57679770e-05, 1.73737888e-05, 1.90444589e-05, 1.95464853e-05,
            1.96534246e-05, 2.01932476e-05, 2.11411082e-05, 2.23337305e-05,2.23086067e-05, 2.15692808e-05,
            2.03738466e-05, 1.86966020e-05,1.77262511e-05, 1.66076057e-05, 1.56260175e-05, 1.49018889e-05,
            1.40255072e-05, 1.36508092e-05, 1.29322568e-05, 1.21598923e-05,1.16423180e-05, 1.10308910e-05,
            1.06868889e-05, 1.00614393e-05,9.37936198e-06, 8.91804495e-06, 8.51265689e-06, 8.39893645e-06,
            8.19183918e-06, 7.98025714e-06, 7.87300425e-06, 7.69113539e-06,7.54182133e-06, 7.23826933e-06,
            6.76030395e-06, 6.00486097e-06,6.96466431e-06, 9.21183232e-06, 7.61517318e-06, 6.18463214e-06,
            6.28987583e-06, 4.98552335e-06, 4.92512830e-06, 4.42123858e-06,3.63720401e-06, 3.41846086e-06,
            2.59899532e-06, 2.14558929e-06,1.76456045e-06, 1.35914120e-06, 1.16952604e-06, 1.40658209e-06};



  for(int iB = 0; iB < nBins; iB++){
    chi2 += pow(Observed[iB] - Ratio[iB],2) / (pow(StatsError[iB],2) + Cov_diag[iB]);
  }
*/
  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "NEOS Chi2: " << chi2 << std::endl;

  return chi2;
}
