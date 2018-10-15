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

  // prompt energy
  Energy = {          1.05554, 1.15959, 1.25667, 1.35333, 1.45272, 1.55443, 1.65394, 1.75384, 1.85395,
                      1.95669, 2.05691, 2.15668, 2.25656, 2.35337, 2.45401, 2.55558, 2.65371, 2.75365,
                      2.85314, 2.95504, 3.05502, 3.15733, 3.25636, 3.35703, 3.45452, 3.55391, 3.65320,
                      3.75569, 3.85493, 3.95334, 4.05519, 4.15477, 4.25459, 4.35448, 4.45603, 4.55511,
                      4.65367, 4.75292, 4.85386, 4.95462, 5.05290, 5.15490, 5.25578, 5.35309, 5.45342,
                      5.55353, 5.65251, 5.75149, 5.85396, 5.95279, 6.05050, 6.15142, 6.25248, 6.35596,
                      6.45501, 6.55417, 6.65510, 6.75543, 6.85126, 6.95296 };
  // neos/dayabay data
  Observed = {        0.97774, 1.03994, 1.03064, 0.96082, 0.95823, 0.96296, 0.96341, 0.97332, 0.98857,
                      0.96860, 0.98643, 0.99329, 1.00274, 0.98689, 1.01509, 1.01616, 0.98262, 0.99360,
                      0.99360, 1.00290, 1.01494, 1.03430, 1.02287, 1.00198, 1.00259, 1.00000, 0.99512,
                      1.01890, 1.01280, 0.98598, 0.99405, 0.99619, 1.00442, 1.01402, 1.01479, 1.00473,
                      1.03171, 1.02576, 1.01128, 1.04253, 1.01280, 1.02454, 1.00869, 1.00488, 0.97561,
                      0.99070, 0.97820, 0.96540, 1.03887, 1.02256, 0.97866, 1.01387, 1.00244, 1.05046,
                      1.03979, 1.03155, 1.01692, 1.03750, 0.99741, 1.05168};
  // neos/dayabay statistical error
  StatsError = {      0.0215, 0.01966, 0.01799, 0.01647, 0.01418, 0.01402, 0.01327, 0.0122, 0.01189,
                      0.01189, 0.01144, 0.01113, 0.01113, 0.01098, 0.01098, 0.01067, 0.01098, 0.01021,
                      0.01082, 0.01051, 0.01067, 0.01097, 0.01128, 0.01128, 0.01159, 0.01174, 0.01174,
                      0.0122, 0.01281, 0.01295, 0.01311, 0.01326, 0.01403, 0.01433, 0.01478, 0.01463,
                      0.01631, 0.01647, 0.01723, 0.01738, 0.01845, 0.01936, 0.01966, 0.02042, 0.02165,
                      0.02332, 0.0247, 0.02591, 0.02775, 0.02912, 0.03125, 0.03339, 0.03445, 0.03902,
                      0.04192, 0.04375, 0.05305, 0.06052, 0.07134, 0.08841};

  // Fluxes
  // NOTE to account for high and  low energy, I incorrectly assume linearity before first point and asymptote to zero for the end
  double fx_energy[] = {1.8, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75,
                        5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 12};
  double fx_u235[] = {  1.48, 1.32, 1.12, 9.15e-1, 7.7e-1, 6.51e-1, 5.53e-1,
                        4.54e-1, 3.643e-1, 2.94e-1, 2.3e-1, 1.79e-1, 1.38e-1,
                        1.1e-1, 8.64e-2, 6.46e-2, 5.1e-2, 3.89e-2, 2.87e-2,
                        2.17e-2, 1.61e-2, 1.14e-2, 7.17e-3, 4.64e-3, 2.97e-3,
                        1.62e-3,0};
  double fx_u238[] = {  1.2724, 1.15, 9.97e-1, 8.55e-1, 7.27e-1, 6.11e-1, 5.06e-1,
                        4.15e-1, 3.36e-1, 2.67e-1, 2.10e-1,  1.63e-1, 1.27e-1,
                        9.69e-2,  7.33e-2, 5.52e-2, 4.14e-2, 3.10e-2, 2.30e-2,
                        1.66e-2,  1.16e-2, 7.85e-3, 5.23e-3, 3.44e-3, 2.19e-3,
                        1.38e-3,0};
  double fx_pu239[] = { 1.208, 1.08, 9.2e-1, 7.19e-1, 6.2e-1, 5.15e-1, 3.98e-1,
                        3.29e-1, 2.61e-1, 1.95e-1, 1.57e-1, 1.13e-1, 8.33e-2,
                        6.13e-2, 4.83e-2, 3.54e-2, 2.92e-2, 1.92e-2, 1.28e-2,
                        9.98e-3, 7.54e-3, 4.98e-3, 3.26e-3, 1.95e-3, 8.47e-4,
                        5.87e-4,0};
  double fx_pu241[] = {1.404, 1.26, 1.08, 8.94e-1, 7.77e-1, 6.41e-1, 5.36e-1,
                        4.39e-1, 3.46e-1, 2.82e-1, 2.2e-1, 1.66e-1, 1.25e-1,
                        9.74e-2, 7.74e-2, 5.58e-2, 4.11e-2, 3.05e-2, 1.98e-2,
                        1.54e-2, 1.09e-2, 7.75e-3, 4.47e-3, 2.9e-3, 1.78e-3,
                        1.06e-3,0};
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
    XSec[i] = dif1.Eval(en+.8);//* 1e-42;
  }

  // Isotope fract (array in order u235 u238 pu239 pu241) for each detector
  DB_f_iso = {{ { .564, .564, .557, .552},
                { .076, .076, .076, .076},
                { .303, .303, .312, .315},
                { .056, .056, .055, .057}}};
  NEOS_f_iso = { 0.655, 0.072, 0.235, 0.038};

  // Efficiency
  double eff_mu[] = { .8255, .8221, .8573, .8571};
  double eff_mult[] = {.9744, .9747, .9757, .9757};
  DB_eff = {eff_mu[0] * eff_mult[0], eff_mu[1] * eff_mult[1], eff_mu[2] * eff_mult[2], eff_mu[3] * eff_mult[3]};

  // Length  (array for each detector in order AD1 AD2 AD3 AD8) for each reactor
  DB_l_d = {{ {362.38, 371.76, 903.47, 817.16, 1353.62, 1265.32},
              {357.94, 368.41, 903.35, 816.90, 1354.23, 1265.89},
              {1332.48, 1358.15, 467.57, 489.58, 557.58, 499.2},
              {1337.43, 1362.88, 472.97, 495.35, 558.71, 501.07}}};

  DB_mass = { 907.f, 916.f, 915.f,  950.f };

  // Smearing matrix
  file.open(dataLoc+"DayaBay_DetectorModel_Etrue_to_Erec.txt");
  for(short i = 0; i < 240; i++){
    for(short j = 0; j < 240; j++){
  	  file >> DB_smearingmatrix[i][j];
    }
  }
  file.close();

  // Generate the smearing matrix for NEOS_H
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

  //
  // Daya Bay 3v
  //
  std::array < std::array < double, 200 >, 4 > DB_Predicted_3n_Unsmeared, DB_Predicted_3n_Smeared;
  for(int iDet = 0; iDet < 4; iDet ++){ // Loop through detectors
    for(int iEn = 0; iEn < 200; iEn++){
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

  dof = nBins;

  // Iniitalize our output tree
  chi2Nt = new OutTree("NEOS");

	if(debug) std::cout << "NEOS initialized. Bins: " << nBins << std::endl;

  return dof;
}

float NEOS::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

	oscContribution oscCon = getOscContributionsNueDis(model);

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
            prob = 1.;
            for(int iContribution = 0; iContribution < 6; iContribution++){
              prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * DB_l_d[iDet][iReac] / (en+.8)),2);
            }
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
          prob = 1.;
          for(int iContribution = 0; iContribution < 6; iContribution++){
            prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * len / (en+.8)),2);
          }
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
  // fucking wooof what a journey.
  std::array < double, nBins > Ratio;
  for(int iB = 0; iB < nBins; iB++){
    Ratio[iB] = NEOS_Ratio[iB]/DB_Ratio[iB];
    //std::cout << "RATIO: " << Ratio[iB] << std::endl;

    // Add stats error to cov Matrix
    //Cov[iB][iB] += pow(StatsError[iB],2);
  }

  // Calc chi2
/*
  for(int iB = 0; iB < nBins; iB++){
    for(int jB = 0; jB < nBins; jB++){
      chi2 += (Observed[iB] - Ratio[iB]) * Cov[iB][jB] * (Observed[jB] - Ratio[jB]);
//      if(chi2 < 0)
//        std:cout << "NEG: " << (Observed[iB] - Ratio[iB]) << " " <<  (Observed[jB] - Ratio[jB]) << " COV: " << Cov[iB][jB] << std::endl;
    }
  }
*/
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

  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "NEOS Chi2: " << chi2 << std::endl;

  return chi2;
}
