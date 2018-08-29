#include "Gallium.h"

int Gallium::Init(std::string dataLoc, Oscillator osc, bool debug){

  // Gallium's nice and easy up here.
  double temp0[] = {1.,0.81,0.95,0.79};
  double temp1[] = {0.1,0.1,0.12,0.1};
  double temp2[] = {0.811,0.813};
  double temp3[] = {0.902,0.098};
  double temp4[] = {70.1,70.3};
  double temp5[] = {0.747,0.742,0.427,0.432};
  double temp6[] = {0.8163,0.0849,0.0895,0.0093};
  double temp7[] = {60.8,61.5,26.7,27.1};

  for(int i = 0; i < 4; i++){
  	obsRatioGal[i] = temp0[i];
  	errorGal[i] = temp1[i];
  	crLinesE[i] = temp5[i];
  	crLinesBr[i] = temp6[i];     // neutrino fraction
  	crLinesXSec[i] = temp7[i];   // neutrino capture xsec
  }
  for(int i = 0; i < 2; i++){
  	arLinesE[i] = temp2[i];
  	arLinesBr[i] = temp3[i];     // neutrino fraction
  	arLinesXSec[i] = temp4[i];   // neutrino capture xsec
  }

	// Here's a little experiment.
	// so, L will go from 0 to the corner. In this way, we'll magically reduce two parameters to only one!
	// we need three L's and their corresponding volume integrals
	double nLength, nInt;
	nLength = 2000; nInt = 600;

	volInt.resize(3, std::vector<double>(nLength));
	for(int i = 0; i < nLength; i++){
		volInt[0][i] = 0.;	volInt[1][i] = 0.;	volInt[2][i] = 0.;
	}

	for(int iex = 0; iex < 3; iex++){
		double maxht = max(Height[iex]-SourceHeight[iex],SourceHeight[iex]);
		double maxlen = sqrt(pow(Radius[iex],2) + pow(maxht,2));

		double ht, dht, rad, dr;
		for(int iHt = 0; iHt < nInt; iHt++){
			ht = (Height[iex] / nInt) * (iHt + .5); dht = Height[iex] / nInt;

			for(int iRd = 0; iRd < nInt; iRd++){
				rad = (Radius[iex] / nInt) * iRd; dr = Radius[iex] / nInt;

				int _ind = floor(sqrt(pow(rad,2) + pow(ht-SourceHeight[iex],2))/(maxlen/float(nLength)));
				volInt[iex][_ind] += dr * dht * 2*TMath::Pi()*rad;
			}
		}
	}

  dof = 4;

  // Iniitalize our output tree
  chi2Nt = new OutTree("Gallium");

	if(debug) std::cout << "Gal initialized. Bins: " << 4 << std::endl;

  return dof;
}

float Gallium::Chi2(Oscillator osc, neutrinoModel model,bool debug){

  float chi2 = 0.f;

  int nLinesCr = 4; int nLinesAr = 2;
  double denominator[4] = {1674.,1675.,580.8,708.4};
  double numerator[4];

	oscContribution oscCon = getOscContributionsNueDis(model);

  for(int iG = 0; iG < nPoints; iG++){
    numerator[iG] = 0;
    if(iG < 3){
		  for(int iCr = 0; iCr < 4; iCr++){
        double maxht = max(Height[iG]-SourceHeight[iG],SourceHeight[iG]);
			  double maxlen = sqrt(pow(Radius[iG],2) + pow(maxht,2));

				for(int iLen = 0; iLen < 2000; iLen++){

					double prob = 1.;
					for(int iContribution = 0; iContribution < 6; iContribution++){
						prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * (iLen+1)/float(2000)*maxlen / crLinesE[iCr]),2);
					}
					numerator[iG] += volInt[iG][iLen] * (1./pow((iLen+1)/float(2000)*maxlen,2)) * prob * crLinesBr[iCr] * crLinesXSec[iCr];
				}
			}
    }
    else if(iG == 3){
      for(int iAr = 0; iAr < 2; iAr++){
				double maxht = max(Height[2]-SourceHeight[2],SourceHeight[2]);
				double maxlen = sqrt(pow(Radius[2],2) + pow(maxht,2));

        for(int iLen = 0; iLen < 2000; iLen++){

					double prob = 1.;
					for(int iContribution = 0; iContribution < 6; iContribution++){
						prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * (iLen+1)/float(2000)*maxlen / arLinesE[iAr]),2);
					}
					numerator[iG] += volInt[2][iLen] * (1./pow((iLen+1)/float(2000)*maxlen,2)) * prob * arLinesBr[iAr] * arLinesXSec[iAr];
				}
      }
    }
	}

  FullData.resize(nPoints);
  Prediction.resize(nPoints);

  for(int iG = 0; iG < nPoints; iG++){
    Prediction[iG] = numerator[iG] / denominator[iG];
    FullData[iG] = obsRatioGal[iG];
  }

  for(int iG = 0; iG < nPoints; iG++){
    chi2 += pow(FullData[iG] - Prediction[iG],2) / pow(errorGal[iG],2);
  }

  // Fill output Tree
  chi2Nt->Fill(chi2, dof, model);

  if(debug)
    std::cout << "Gallium Chi2: " << chi2 << std::endl;

  return chi2;
}
