#include "DANSS.h"

/*
Davio Cianci Aug 26th

Fit method: https://arxiv.org/pdf/1709.04294.pdf
DANSS paper: https://arxiv.org/pdf/1804.04046.pdf
Energy Resolution: https://arxiv.org/pdf/1412.0817.pdf

*/


int DANSS::Init(std::string dataLoc, Oscillator osc, bool debug){

  Observed.resize(nBins); // observed event ratio of down/up detectors
  Energy.resize(nBins);   // Positron energy;
  StatsError.resize(nBins);

  double temp1[] = {0.7200, 0.7219, 0.7306, 0.7148, 0.7132, 0.7094, 0.7123, 0.7110, 0.6911, 0.7011, 0.7158, 0.6996,
                  0.7054, 0.6904, 0.6880, 0.7157, 0.7168, 0.7422, 0.7272, 0.7616, 0.7538, 0.6994, 0.6899, 0.6694};
  double temp2[] = { 0.0127662, 0.0103213, 0.00884577, 0.00879326, 0.00846355, 0.00879573, 0.00850462, 0.00895212,
                  0.0099994, 0.00960843, 0.00890215, 0.0109418, 0.0109048, 0.0125904, 0.0138302, 0.0132009, 0.0141713,
                  0.0148419, 0.0182969, 0.0192751, 0.0235301, 0.0299239, 0.0361098, 0.043911};

  for(int i = 0; i < nBins; i++){
    Observed[i] = temp1[i];
    StatsError[i] = temp2[i];
  }

  EnergySmearing.resize(48,std::vector<double>(48));
  // Load in smearing matrix
  ifstream file;
	file.open(dataLoc+"danss_smearingmatrix.txt");
	for(short i = 0; i < 48; i++)
    for(short j = 0; j < 48; j++)
		  file >> EnergySmearing[i][j];
	file.close();

  dof = nBins;

  // Iniitalize our output tree
  chi2Nt = new OutTree("DANSS");

	if(debug) std::cout << "DANSS initialized. Bins: " << nBins << std::endl;

  return dof;
}

float DANSS::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

	oscContribution oscCon = getOscContributionsNueDis(model);

  std::vector < double > Numerator, Denominator, ProbRatio, RatioSmeared;
  Numerator.resize(48);        // Down detector prob
  Denominator.resize(48);      // Up detector prob
  ProbRatio.resize(48);
  RatioSmeared.resize(48);

  // Calculate ratio for positron spectrum and then apply gaussian smearing
  int binAvg = 50;
  int lenAvg = 10;
  for(int ei = 0; ei < 48; ei++){
    Numerator[ei] = 0;
    Denominator[ei] = 0;
    // Average across each bin to catch fast oscillations (particularly in low bins)
    for(int i = 0; i < binAvg; i++){
      double en = (.25 * ei) + i / float(binAvg) * .25;
      // Also, average across the length of the detector!
      for(int i = 0; i < lenAvg; i++){
        double len_down = osc.RanGen.Gaus(L0_down,4./3.);
        double len_up = osc.RanGen.Gaus(L0_up,4./3.);
        double probDown = 1.;
        double probUp = 1.;
        for(int iContribution = 0; iContribution < 6; iContribution++){
          probDown += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * len_down / (en+1.8)),2);
          probUp += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * len_up / (en+1.8)),2);
        }
        Numerator[ei] += probDown / pow(12.7,2);
        Denominator[ei] += probUp / pow(10.7,2);
      }
    }
  }

  // Now, we smear everything into proper bins!
  for(int ei = 0; ei < 48; ei++){
    RatioSmeared[ei] = 0;
    for(int j = 0; j < 48; j++){
      RatioSmeared[ei] += Numerator[j]/Denominator[j] * EnergySmearing[ei][j];
    }
  }

  // Lastly, just take out the bit that we're interested in, from positron energy 1-77MeV
  for(int iB = 0; iB < 24; iB++)
    ProbRatio[iB] = RatioSmeared[iB+4];

  double systerror = .02;
  for(int iB = 0; iB < nBins; iB++){
    chi2 += pow((Observed[iB] - ProbRatio[iB]) / sqrt(pow(StatsError[iB],2) + pow(systerror*ProbRatio[iB],2)),2);
    //std::cout << "A: " << ProbRatio[iB] <<  std::endl;
  }

  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "DANSS Chi2: " << chi2 << std::endl;

  return chi2;
}
