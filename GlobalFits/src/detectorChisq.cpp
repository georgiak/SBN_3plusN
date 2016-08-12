/* ------------------------------------------//
Created by Davio Cianci and Georgia Karagiorgi
Jan 25th, 2016

------------------------------------------// */
#include "TStopwatch.h"


TMatrixT <double> cov;      // inverted cov matrix
TMatrixT <double> covMatrix;
TMatrixT <double> full_covMatrix;
double cov_det;
std::vector <double> _signal;
std::vector <double> _fullData;
std::vector <double> _prediction;
std::vector <double> _totalError;
TMinuit *gMinuit;
minPack myMin;

void myMinInit(){
	gMinuit = new TMinuit(6);
}

chisqStruct getChi2Boone(neutrinoModel model, boonePackage pack, bool nubar){

    int nBins = 11;
    int nBins_mu = 8;
    int nFOscEvts = pack.nFOscEvts;

    // Initialize the result!
    chisqStruct result;
    result.zero();

    // Initialize contributions from the oscillation probability
    oscContribution oscCont;

    _signal.resize(nBins);
    _fullData.resize(nBins + nBins_mu);
    _prediction.resize(nBins + nBins_mu);

    double minEBins[nBins];
    double maxEBins[nBins];
    double ETru[nFOscEvts], LTru[nFOscEvts];

    // Initialize vars
    for(int iB = 0; iB < nBins; iB++){
        _signal[iB] = 0.;
        _fullData[iB] = pack.NueData[iB];
        _prediction[iB] = pack.NueBgr[iB];
    }
    for(int iB = 0; iB < nBins_mu; iB++){
        _fullData[iB + nBins] = pack.NumuData[iB];
        _prediction[iB + nBins] = pack.Numu[iB];
    }
    full_covMatrix.ResizeTo(nBins + nBins + nBins_mu, nBins + nBins + nBins_mu);
    full_covMatrix.Zero();
    covMatrix.ResizeTo(nBins + nBins_mu, nBins + nBins_mu);
    covMatrix.Zero();

	// Get oscillation probability contributions
	oscCont = getOscContributionsNueApp(model, nubar, true);

	for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
            minEBins[iB] = pack.EnuQE[iB];
            maxEBins[iB] = pack.EnuQE[iB+1];

            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
                // Get prediction signal by multiplying osc prob by weight of each event

                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

                for(int iContribution = 0; iContribution < 6; iContribution++){

                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
                }
            }
        }
    }

    // Divide signal prediction by the number of fullosc events
    for(int iB = 0; iB < nBins; iB++){
        _signal[iB] /= double(nFOscEvts);
    }

    // Now, scale the fractional cov matrix to our signal and prediction vectors
    for(int iB = 0; iB < nBins + nBins + nBins_mu; iB++){
        for(int jB = 0; jB < nBins + nBins + nBins_mu; jB++){
            if(iB < nBins && jB < nBins){
                full_covMatrix(iB,jB) = pack.full_fractCovMatrix[iB][jB]*_signal[iB]*_signal[jB];
                // Add Stat error of signal prediction
                if(iB == jB){
                    full_covMatrix(iB,jB) += _signal[iB];
                }
            }
            else if(iB < nBins && jB >= nBins){
                full_covMatrix(iB,jB) = pack.full_fractCovMatrix[iB][jB]*_signal[iB]*_prediction[jB-nBins];
            }
            else if(iB >= nBins && jB < nBins){
                full_covMatrix(iB,jB)= pack.full_fractCovMatrix[iB][jB]*_prediction[iB-nBins]*_signal[jB];
            }
            else if(iB >= nBins && jB >= nBins){
                full_covMatrix(iB,jB) = pack.full_fractCovMatrix[iB][jB]*_prediction[iB-nBins]*_prediction[jB-nBins];
            }
        }
    }

    // Now, collapse our 3x3 matrix to a 2x2
    for(int iB = 0; iB < nBins + nBins_mu; iB++){
        for(int jB = 0; jB < nBins + nBins_mu; jB++){
            if(iB < nBins && jB < nBins){
                covMatrix(iB,jB) = full_covMatrix[iB][jB] + full_covMatrix[iB + nBins][jB] + full_covMatrix[iB][jB + nBins] + full_covMatrix[iB + nBins][jB + nBins];
            }
            else if(iB < nBins && jB >= nBins){
                covMatrix(iB,jB) = full_covMatrix[iB][jB + nBins] + full_covMatrix[iB + nBins][jB + nBins];
            }
            else if(iB >= nBins && jB < nBins){
                covMatrix(iB,jB) = full_covMatrix[iB + nBins][jB] + full_covMatrix[iB + nBins][jB + nBins];
            }
            else if(iB >= nBins && jB >= nBins){
                covMatrix(iB,jB) = full_covMatrix[iB + nBins][jB + nBins];
            }
        }
    }

    // Now, let's invert the covariance matrix
    cov.ResizeTo(nBins + nBins_mu, nBins + nBins_mu);
    cov = covMatrix.Invert();

    for(int iB = 0; iB < nBins; iB++){
        _prediction[iB] += _signal[iB];
    }

    // Now, let's calculate the determinant of the cov matrix
    cov_det = cov.Determinant();

    // Finally, let's put everything together and calculate the chisq
    for(int iB = 0; iB < nBins + nBins_mu; iB++){
        for(int jB = 0; jB < nBins + nBins_mu; jB++){
            result.chi2 += (_fullData[iB]-_prediction[iB])*cov(iB,jB)*(_fullData[jB]-_prediction[jB]);
        }
    }

    result.chi2_det = result.chi2_det + cov_det;
    return result;
}
// Gallium
chisqStruct getChi2Gallium(neutrinoModel model, galPackage pack){

	chisqStruct result;
    result.zero();

    const int nPoints = 4; int nLinesCr = 4; int nLinesAr = 2;

    double radiusGallex = 1.9;  double heightGallex = 5.0;  double sourceHeightGallex[2] = {2.7,2.38};
    double radiusSage = 0.7;    double heightSage = 1.47;   double sourceHeightSage = 0.72;
    double denominator[4] = {1674.,1675.,580.8,708.4};
    double numerator[4];

	double Radius[3] = {radiusGallex, radiusGallex, radiusSage};
	double Height[3] = {heightGallex, heightGallex, heightSage};
	double SourceHeight[3] = {sourceHeightGallex[0], sourceHeightGallex[1], sourceHeightSage};

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
						prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * (iLen+1)/float(2000)*maxlen / pack.crLinesE[iCr]),2);
					}
					numerator[iG] += pack.volInt[iG][iLen] * (1./pow((iLen+1)/float(2000)*maxlen,2)) * prob * pack.crLinesBr[iCr] * pack.crLinesXSec[iCr];
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
						prob += oscCon.aEE[iContribution] * pow(sin(1.267 * oscCon.dm2[iContribution] * (iLen+1)/float(2000)*maxlen / pack.arLinesE[iAr]),2);
					}
					numerator[iG] += pack.volInt[2][iLen] * (1./pow((iLen+1)/float(2000)*maxlen,2)) * prob * pack.arLinesBr[iAr] * pack.arLinesXSec[iAr];
				}
        	}
    	}
	}

    _fullData.resize(nPoints);
    _prediction.resize(nPoints);

    for(int iG = 0; iG < nPoints; iG++){
        _prediction[iG] = numerator[iG] / denominator[iG];
        _fullData[iG] = pack.obsRatioGal[iG];
    }

    for(int iG = 0; iG < nPoints; iG++){
        result.chi2 += pow(_fullData[iG] - _prediction[iG],2) / pow(pack.errorGal[iG],2);
    }

	std::cout << result.chi2 << std::endl;

    return result;
}
// Atmospheric - though this barely counts
chisqStruct getChi2Atm(neutrinoModel model, atmPackage pack){

    chisqStruct result;
    result.zero();

    int dmuVecMaxDim = 101;
    double dmu = 0;

    ROOT::Math::Interpolator dif(dmuVecMaxDim);
    dif.SetData(dmuVecMaxDim,pack.dmuVec,pack.dchi2Vec);

    // First, let's find dmu such that dmu**2 - dmu + A = 0
    double A = (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2)) + pow(model.Um[0]*model.Um[1],2) + pow(model.Um[0]*model.Um[2],2) + pow(model.Um[1]*model.Um[2],2);
    dmu = (1 - sqrt(1. - 4*A))/2.;

    // For this, everything has been figured out, so we just interpolate an array of chi2's to get our result as a function of dmu
	//std::cout << pack.dmuVec[0] << " " << dmu << " " << pack.dmuVec[dm2VecMaxDim-1] << std::endl;
	result.chi2 = dif.Eval(dmu);

    return result;
}
// NUMI stuff
chisqStruct getChi2Numi(neutrinoModel model, numiPackage pack){

    int nBins = 10;
    int nFOscEvts = 3323;

    // Initialize the result!
    chisqStruct result;
    result.zero();

    // Initialize contributions from the oscillation probability
    oscContribution oscCont;

    _signal.resize(nBins);
    _totalError.resize(nBins);

    double minEBins[nBins];
    double maxEBins[nBins];
    double ETru[nFOscEvts], LTru[nFOscEvts];

    for(int i = 0; i < nBins; i++){
        _signal[i] = 0.;
        _totalError[i] = 0;
    }

	oscCont = getOscContributionsNueApp(model, false, true);

    // Loop over FOsc Events
    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){
        // Loop over energy bins
        for(int iN = 0; iN < nBins; iN++){
            minEBins[iN] = pack.EnuQE[iN];
            maxEBins[iN] = pack.EnuQE[iN+1];

            //std::cout << minEBins[iN] << " // " << pack.EnuQE[iN] << " " << pack.FOsc_EnuQE[iFOsc] << " " << pack.EnuQE[iN+1] << " // " << maxEBins[iN] << std::endl;
            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iN] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iN]){
                // Now that we're in the right energy range, get the predicted signal
                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

                for(int iContribution = 0; iContribution < 6; iContribution++){
                    //std::cout << pack.FOsc_weight[iFOsc] << " " << oscCont.aMuE[iContribution] << " " << oscCont.dm2[iContribution] << " " << LTru[iFOsc] << " " << ETru[iFOsc] << std::endl;
                    _signal[iN] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow((sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]/ETru[iFOsc])),2)
                        + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2.*oscCont.dm2[iContribution]*LTru[iFOsc]/ETru[iFOsc]);
                }
            }
        }
    }

    // Now fill up that error
    for(int iN = 0; iN < nBins; iN++){
        _totalError[iN] = pow(pack.NueBgr_error[iN],2) + pow(pack.FOsc_fracError[iN]*_signal[iN],2) + _signal[iN] + pack.NueBgr[iN];
	}

    // Now actually calculate the chisq
    for(int iN = 0; iN < nBins; iN++){
        result.chi2 += pow((pack.NueData[iN] - pack.NueBgr[iN] - _signal[iN]),2) / _totalError[iN];
    }

    return result;
}
// MiniBoone disappearance
chisqStruct getChi2MBDis(neutrinoModel model, booneDisPackage pack){

    const int nBins = 16;

    // Initialize the result!
    chisqStruct result;
    result.zero();

	oscContribution oscCont;

    double minEBins[nBins], maxEBins[nBins];
    double dtIntegral = 0.;     double MCIntegral = 0.;
	double ETru, LTru;
	double FOsc_EnuQE, FOsc_EnuTrue, FOsc_LnuTrue, FOsc_weight;

    _signal.resize(nBins);
    _prediction.resize(nBins);

    for(int iB = 0; iB < nBins; iB ++){
        _signal[iB] = 0;
        _prediction[iB] = 0;
        dtIntegral += pack.NumuData[iB];
    }

    covMatrix.ResizeTo(nBins, nBins);
    covMatrix.Zero();

	oscCont = getOscContributionsNumuDis(model);

	for(int iEvt = 0; iEvt < pack.nFOscEvts; iEvt++){
		//if(iEvt%5==0) continue;

		for(int iB = 0; iB < nBins; iB++){
        	minEBins[iB] = pack.EnuQE[iB];
        	maxEBins[iB] = pack.EnuQE[iB+1];

        	if(pack.FOsc_EnuQE[iEvt] > minEBins[iB] && pack.FOsc_EnuQE[iEvt] < maxEBins[iB]){
            	// Get predicted signal by multiplying the osc prob by the weight of each event!
            	ETru = pack.FOsc_EnuTrue[iEvt];
            	LTru = pack.FOsc_LnuTrue[iEvt];

            	// No-osc signal prediction
        		_prediction[iB] += pack.FOsc_weight[iEvt];

        		for(int iContribution = 0; iContribution < 6; iContribution++){
            		_signal[iB] += pack.FOsc_weight[iEvt]*oscCont.aMuMu[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru / ETru),2);
            	}
			}
		}
	}

    // Normalize signal prediction to data
    for(int iB = 0; iB < nBins; iB++){
        MCIntegral += (_prediction[iB] + _signal[iB]);
    }
    for(int iB = 0; iB < nBins; iB++){
        _prediction[iB] = (_prediction[iB] + _signal[iB]) * dtIntegral/MCIntegral;
    }

    for(int iB = 0; iB < nBins; iB++){
        for(int jB = 0; jB < nBins; jB++){
            covMatrix[iB][jB] = pack.full_fractCovMatrix[iB][jB] * _prediction[iB] * _prediction[jB];
            // Add statistical error of signal prediction
            if(iB == jB){
                covMatrix[iB][jB] += _prediction[iB];
            }
        }
    }

    // Now, let's invert this bad boy
    cov.ResizeTo(nBins,nBins);
    cov = covMatrix.Invert();

    // Finally, let's put everything together and calculate the chisq
    for(int iB = 0; iB < nBins; iB++){
        for(int jB = 0; jB < nBins; jB++){
            result.chi2 += (pack.NumuData[iB] - _prediction[iB]) * cov[iB][jB] * (pack.NumuData[jB] - _prediction[jB]);
        }
    }

    return result;
}
// Minos CC
double getChi2MinosSign(neutrinoModel model, const int nBins, double * EnuQE, double * NumubarData, double * NumubarBkg, double * fracError, double * dataErr){

    double chi2 = 0.;

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
    _signal.resize(nBins);
    _prediction.resize(nBins);
    _totalError.resize(nBins);

    for(int i = 0; i < nBins; i++){
        _signal[i] = 0;
        _prediction[i] = 0;
        _totalError[i] = 0;
    }

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

        _prediction[iM] = NumubarBkg[iM] * (1 - 4 * ((1 - pow(Um3,2) - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(Um3,2) * sin31Sq + pow(model.Um[0],2) * sin41Sq + pow(model.Um[1],2) * sin51Sq + pow(model.Um[2],2) * sin61Sq)
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

        _signal[iM] = (1 - 4* ((1 - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) * sin41Sq + pow(model.Um[1],2) * sin51Sq + pow(model.Um[2],2) * sin61Sq)
                + pow(model.Um[1],2) * pow(model.Um[0],2) * sin54Sq + pow(model.Um[2],2) * pow(model.Um[1],2) * sin65Sq + pow(model.Um[2],2) * pow(model.Um[0],2) * sin64Sq));

        _prediction[iM] /= _signal[iM];

        if(_prediction[iM] < 0)     _prediction[iM] = 0;

        // Now, fill up the total error while we're still loopin' around
        double error_from_atm_mixing = NumubarBkg[iM] * sin22Theta_atm_err * pow(sin(1.27 * dm2Atm * LTrue / ETrue[iM]),2);
        _totalError[iM] = pow(dataErr[iM],2) + pow(fracError[iM] * _prediction[iM],2) + pow(error_from_atm_mixing,2);
	}

    // Now, calculate the right-sign chi2
    for(int iM = 0; iM < nBins; iM++){
        chi2 += (NumubarData[iM] - _prediction[iM]) * (NumubarData[iM] - _prediction[iM]) / _totalError[iM];
	}

    return chi2;
}
chisqStruct getChi2Minos(neutrinoModel model, minosPackage pack){
    chisqStruct result;
    result.zero();

    // MINOS is super annoying, so it's split up into two parts:
    // Right Sign
    double chi2Temp1 = getChi2MinosSign(model, 12, pack.EnuQE, pack.NumubarData, pack.NumubarBkg, pack.fracError, pack.dataErr);
	// Wrong Sign
    double chi2Temp2 = getChi2MinosSign(model, 13, pack.EnuQE_ws, pack.NumubarData_ws, pack.NumubarBkg_ws, pack.fracError_ws, pack.dataErr_ws);

    result.chi2 = chi2Temp1 + chi2Temp2;
    return result;
}
// Karmen and LSND both use this guy! (still a bit of weirdness, maybe. KARMEN always gives very small signal, but maybe the detector is just not super sensitive. Look into it later)
chisqStruct getLogLikelihood(neutrinoModel model, int nBins, sinSqPackage pack){

    chisqStruct result;
    result.zero();
    oscContribution oscCont;

    ROOT::Math::Interpolator dif(dm2VecMaxDim);
    double sinSq, sinSq2;
    double sinSqDeltaVec[dm2VecMaxDim], sinSqDeltaVec2[dm2VecMaxDim];
    double lt1, lt2, lt3, pred;
	double lsndDm2Vec[601];

     _signal.resize(nBins);
    for(int i = 0; i < nBins; i++) _signal[i] = 0.;

	oscCont = getOscContributionsNueApp(model,true,true);
    for(int iL = 0; iL < nBins; iL ++){
        for(int k = 0; k < 601; k ++){
            sinSqDeltaVec[k] = pack.sinSqDeltaGrid[k][iL];
            sinSqDeltaVec2[k] = pack.sinSqDeltaGrid2[k][iL];
			if(pack.dm2Vec.size() == 0) lsndDm2Vec[k] = dm2Vec[k];
			else lsndDm2Vec[k] = pack.dm2Vec[k];
        }

        for(int iContribution = 0; iContribution < 6; iContribution++){
            // Now, add the latest contribution to the predicted signal vector
            dif.SetData(dm2VecMaxDim,lsndDm2Vec,sinSqDeltaVec);
            if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
            else    sinSq = dif.Eval(oscCont.dm2[iContribution]) * pack.norm;
            dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec2);
            if(oscCont.dm2[iContribution] == 0.)    sinSq2 = 0;
            else    sinSq2 = dif.Eval(oscCont.dm2[iContribution]) * pack.norm;

			_signal[iL] += oscCont.aMuE[iContribution] * sinSq + oscCont.aMuE_CPV[iContribution] * sinSq2;
        }
    }

    // Now, using the signal vector, use the log-likelihood method to get the effective chisq
    for(int iL = 0; iL < nBins; iL++){
        pred = _signal[iL] + pack.bkg[iL];
        lt1 = pred - pack.observed[iL];
        if(pred > 0) lt2 = pack.observed[iL] * log(pred);
        else lt2 = 0;
        if(pack.observed[iL] > 0) lt3 = pack.observed[iL] * log(pack.observed[iL]);
        else lt3 = 0;

        result.chi2 += 2 * (lt1 - lt2 + lt3);
    }

    return result;
}
// NOMAD, I suppose
chisqStruct getChi2Nomad(neutrinoModel model, nomadPackage pack){

	TStopwatch *watch = new TStopwatch;
	//CLOCKER
	watch->Stop();	std::cout << "calculate signal" << std::endl;
	watch->Print("u");
	watch->Start();

	chisqStruct result;
    result.zero();
    oscContribution oscCont;

    _signal.resize(30);

    ROOT::Math::Interpolator dif(dm2VecMaxDim);
    const int maxEnergyBins = 10; const int maxRadialBins = 3;
    double l1 = 0.422; double l2 = .837;
    double energyMin[] = {3.,12.,16.,20.,25.,30.,35.,40.,50.,100.};
    double energyMax[] = {12.,16.,20.,25.,30.,35.,40.,50.,100.,170.};

    double numerator, denominator, sinSq, sinSq2;
    double sinSqDeltaVec[dm2VecMaxDim], sinSqDeltaVec2[dm2VecMaxDim];

    // Loop over energy bins
    for(int ii = 0; ii < maxEnergyBins; ii++){
        // Loop over radial bins
        for(int jj = 0; jj < maxRadialBins; jj++){
            int iN = maxEnergyBins*jj + ii;

            for(int k = 0; k < dm2VecMaxDim; k++){
                sinSqDeltaVec[k] = pack.sinSqDeltaGrid[k][iN];
                sinSqDeltaVec2[k] = pack.sinSqDeltaGrid2[k][iN];
            }

            // First, let's do Nue appearance
            numerator = 0.;
            oscCont = getOscContributionsNueApp(model, false, true);

            for(int iContribution = 0; iContribution < 6; iContribution++){
                // Now, add the latest contribution to the predicted signal vector

                dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);
                if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
                else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
                dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec2);
                if(oscCont.dm2[iContribution] == 0.)    sinSq2 = 0;
                else    sinSq2 = dif.Eval(oscCont.dm2[iContribution]);

                numerator += oscCont.aMuE[iContribution]*sinSq + oscCont.aMuE_CPV[iContribution] * sinSq2;
            }
			//std::cout << "Numerator: " << numerator << std::endl;

            // Now, muon disappearance
            denominator = (l2 - l1) * (energyMax[ii] - energyMin[ii]);
            oscCont = getOscContributionsNumuDis(model);

            for(int iContribution = 0; iContribution < 6; iContribution++){
                // Now, add the latest contribution to the predicted signal vector
                dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);
                if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
                else sinSq = dif.Eval(oscCont.dm2[iContribution]);
                denominator += oscCont.aMuMu[iContribution]*sinSq;
            }
            if(denominator < 0) std::cout << "Got a problem with NOMAD. Denominator is negative!" << std::endl;
			//std::cout << "Denominator: " << denominator << std::endl;

			// Electron to muon ratio
            _signal[iN] = pack.norm[iN] * (numerator / denominator);
            //std::cout << pack.norm[iN] << " " << numerator << " " << denominator << " " << _signal[iN] << std::endl;
        }
    }

    for(int i = 0; i < 30; i++){
        for(int j = 0; j < 30; j++){
            if(i == j)
                result.chi2 += (pack.observed[i] - (pack.bkg[i] + _signal[i])) * pack.sigmaRemu[i][j] * (pack.observed[j] - (pack.bkg[j] + _signal[j]));
        }
    }

    return result;
}

chisqStruct getChi2CCFR(neutrinoModel model, ccfrPackage pack){

    chisqStruct result;
    result.zero();
    oscContribution oscCont;

    const int maxEnergyBins = 18;

    ROOT::Math::Interpolator dif(dm2VecMaxDim);
    double sinSqDeltaVec[dm2VecMaxDim];
    double nFront[maxEnergyBins], nBack[maxEnergyBins];
    double nFront_noOsc[maxEnergyBins], nBack_noOsc[maxEnergyBins];
    double sinSq, ratio_theor[maxEnergyBins];
	// Get the oscillation contributions - we've got a disappearance case here, looks like
	oscCont = getOscContributionsNumuDis(model);

    // Front Detector
    for(int iC = 0; iC < maxEnergyBins; iC++){
        for(int k = 0; k < dm2VecMaxDim; k++){
            sinSqDeltaVec[k] = pack.sinSqDeltaGrid_front[k][iC];
        }
        dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);

		nFront[iC] = pack.noOscGrid[iC][0];

        for(int iContribution = 0; iContribution < 6; iContribution++){
			if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
			else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
            nFront[iC] += oscCont.aMuMu[iContribution] * sinSq;
        }
        nFront[iC] *= pack.m_front;
        nFront_noOsc[iC] = pack.m_front * pack.noOscGrid[iC][0];
    }

    // Back Detector
    for(int iC = 0; iC < maxEnergyBins; iC++){
        for(int k = 0; k < dm2VecMaxDim; k++){
            sinSqDeltaVec[k] = pack.sinSqDeltaGrid_back[k][iC];
        }
        dif.SetData(dm2VecMaxDim,dm2Vec,sinSqDeltaVec);

        nBack[iC] = pack.noOscGrid[iC][1];

        for(int iContribution = 0; iContribution < 6; iContribution++){
			if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
			else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
            nBack[iC] += oscCont.aMuMu[iContribution] * sinSq;
        }
        nBack[iC] *= pack.m_back;
        nBack_noOsc[iC] = pack.m_back * pack.noOscGrid[iC][1];
    }

    // Get the predicted ratio
    for(int iC = 0; iC < maxEnergyBins; iC++){
        ratio_theor[iC] = (nBack[iC]/nFront[iC])/(nBack_noOsc[iC]/nFront_noOsc[iC]);
		//std::cout << "ratio theor: " << ratio_theor[iC] << std::endl;
    }

    // Calculate the chisq
    for(int iC = 0; iC < maxEnergyBins; iC++){
        for(int jC = 0; jC < maxEnergyBins; jC++){
            result.chi2 += (pack.observed[iC] - ratio_theor[iC])*pack.sigmaRatio[iC][jC]*(pack.observed[jC] - ratio_theor[jC]);
		}
    }

    return result;
}

chisqStruct getChi2CDHS(neutrinoModel model, cdhsPackage pack){

    chisqStruct result;
    result.zero();
    oscContribution oscCont;

    const int maxEnergyBins = 15;

    ROOT::Math::Interpolator dif(dm2VecMaxDim);
    double sinSqDeltaVec[dm2VecMaxDim];
    double nFront[maxEnergyBins], nBack[maxEnergyBins];
    double nFront_noOsc[maxEnergyBins], nBack_noOsc[maxEnergyBins];
    double sinSq, ratio_theor[maxEnergyBins];
	double cdhsDm2Vec[601];

	// Get the oscillation contributions - we've got a disappearance case here
	oscCont = getOscContributionsNumuDis(model);

    // Front Detector
    for(int iC = 0; iC < maxEnergyBins; iC++){

        for(int k = 0; k < 601; k++){
			cdhsDm2Vec[k] = pack.dm2Vec[k];
            sinSqDeltaVec[k] = pack.sinSqDeltaGrid_front[k][iC];
        }

        dif.SetData(601,cdhsDm2Vec,sinSqDeltaVec);

        nFront[iC] = pack.noOscGrid[iC][0];

        for(int iContribution = 0; iContribution < 6; iContribution++){
			if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
			else    {
				sinSq = dif.Eval(oscCont.dm2[iContribution]);
			}
            nFront[iC] += oscCont.aMuMu[iContribution] * sinSq;
        }
        nFront[iC] *= pack.m_front[iC];
        nFront_noOsc[iC] = pack.m_front[iC] * pack.noOscGrid[iC][0];
    }

    // Back Detector
    for(int iC = 0; iC < maxEnergyBins; iC++){
        for(int k = 0; k < 601; k++){
            sinSqDeltaVec[k] = pack.sinSqDeltaGrid_back[k][iC];
        }
        dif.SetData(601,cdhsDm2Vec,sinSqDeltaVec);

        nBack[iC] = pack.noOscGrid[iC][1];

        for(int iContribution = 0; iContribution < 6; iContribution++){
			if(oscCont.dm2[iContribution] == 0.)    sinSq = 0;
			else    sinSq = dif.Eval(oscCont.dm2[iContribution]);
            nBack[iC] += oscCont.aMuMu[iContribution] * sinSq;
        }
        nBack[iC] *= pack.m_back[iC];
        nBack_noOsc[iC] = pack.m_back[iC] * pack.noOscGrid[iC][1];
    }

    // Get the predicted ratio
    for(int iC = 0; iC < maxEnergyBins; iC++){
        ratio_theor[iC] = (nBack[iC]/nFront[iC])/(nBack_noOsc[iC]/nFront_noOsc[iC]);
    }

    // Calculate the chisq
    for(int iC = 0; iC < maxEnergyBins; iC++){
        for(int jC = 0; jC < maxEnergyBins; jC++){
            result.chi2 += (pack.observed[iC] - ratio_theor[iC])*pack.sigmaRatio[iC][jC]*(pack.observed[jC] - ratio_theor[jC]);
        }
    }

    return result;
}


// These guys below all use minuit, so they're trouble

// Chooz, also more or less working
void fcnChooz(int &npar, double *gin, double &fval, double  *xval, int iflag){

	double *_gin = gin;	int &_npar = npar;	int _iflag = iflag;

	double alpha = xval[0];

	double Enu, prob, pNueNue, chisq;
	double pNueNuevec[14];

	oscContribution oscCon = getOscContributionsNueDis(myMin.model);
	for(int iChooz = 0; iChooz < 14; iChooz ++){
		Enu = myMin.cPack.energy[iChooz] + 1.8;
		prob = 1;

		// At dm2's this high, nuebar survival probability averages out
		pNueNue = .5 * alpha;

		for(int iCon = 0; iCon < 6; iCon++){
			prob += oscCon.aEE[iCon] * pNueNue;
		}
		pNueNuevec[iChooz] = prob;
	}

	chisq = 0;
	for(int i = 0; i < 14; i++){
		for(int j = 0; j < 14; j++){
			chisq += (myMin.cPack.observed[i] - myMin.cPack.noOsc[i] * pNueNuevec[i]) * myMin.cPack.sigmaX[i][j] * (myMin.cPack.observed[j] - myMin.cPack.noOsc[j] * pNueNuevec[j]);
		}
	}
	// Add the normalization uncertainty
	chisq += pow((alpha-1.)/myMin.cPack.sigmaAlpha,2);
	fval = chisq;

	return;
}
chisqStruct getChi2Chooz(neutrinoModel model, choozPackage pack){

	chisqStruct result;
    result.zero();

    double chisq = 0;

	myMin.cPack = pack;
    myMin.model = model;

	double arglis[2];
    int ierflag;
    arglis[0] = -1.;
    gMinuit->SetFCN(fcnChooz);
    // Okay, let's get this minuit garbage started
    gMinuit->mnexcm("SET PRI",arglis,1,ierflag);
	gMinuit->mnexcm("SET NOWARNING",arglis,0,ierflag);
	arglis[0] = 1.;
	gMinuit->mnexcm("SET STR",arglis,1,ierflag);
    // Now, clear parameters and set 'em anew
    gMinuit->mnexcm("CLE",arglis,1,ierflag);
	gMinuit->mnparm(0,TString("alpha"),1.,pack.sigmaAlpha,0.5,1.5,ierflag);
	// With these high dm2's, CHOOZ oscillations average out so we don't need to fit energy scale factor g
	arglis[0] = 5.;
	gMinuit->mnexcm("CALL", arglis,1,ierflag);
	double disttomin, errdef;   int npari, nparx, istat;
    gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
	result.chi2 = chisq;

  	if(chisq > 150)
  		return result;

    // If we're still here, minimize the fucker.
    arglis[0] = 10000.;
    arglis[1] = 100;
    gMinuit->mnexcm("MINI",arglis,2,ierflag);
    gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
	result.chi2 = chisq;

	return result;
}

// lsnd_karmen_xsec thing
double getSinSqTerm(double dm2, double E1, double E2, double l1, double l2){

	double k = 1.27 * dm2;
	double a11 = 2. * k * l1/E1;
	double a12 = 2. * k * l1/E2;
	double a21 = 2. * k * l2/E1;
	double a22 = 2. * k * l2/E2;

	double sinSqTerm = -4 * E2 * k * l2 - 4 * pow(k,2) * sineInt(a21) * pow(l2,2)
			+ 4. * pow(k,2) * sineInt(a22) * pow(l2,2)
			+ 4. * E1 * k * l2 - sin(a21) * pow(E1,2) + sin(a22) * pow(E2,2)
			- 2. * k * cos(a21) * E1 * l2 + 2. * k * cos(a22) * E2 * l2
			+ 4. * E2 * k * l1 + 4. * pow(k,2) * sineInt(a11) * pow(l1,2)
			- 4. * pow(k,2) * sineInt(a12) * pow(l1,2) - 4. * E1 * k * l1 + sin(a11) * pow(E1,2)
			- sin(a12) * pow(E2,2) + 2. * k * cos(a11) * E1 * l1 - 2. * k * cos(a12) * E2 * l1;

	sinSqTerm = -sinSqTerm / (8. * k) / (E2 - E1) / (l2 - l1);

	// Now, check energy resolution and return .5 if necessary
	double EAvg = (E1 + E2)/2.;	double lAvg = (l1 + l2)/2.;
	double n = 1.27 * dm2 * lAvg / (EAvg * TMath::Pi());
	double ENext = 1.27 * dm2 * lAvg / ((n + .5) * TMath::Pi());
	if(n > 2 && abs(EAvg - ENext)/EAvg < myMin.xPack.ESigma)
		sinSqTerm = 0.5;

	return sinSqTerm;
}
double osc_int(double E1, double E2, double l1, double l2){

	myMin.model.difference();

	double sin41sq,sin51sq,sin54sq,sin61sq,sin64sq,sin65sq;

	if(myMin.model.dm41Sq != 0)
		sin41sq = getSinSqTerm(myMin.model.dm41Sq,E1,E2,l1,l2);
	else
		sin41sq = 0;

	if(myMin.model.dm51Sq != 0){
		sin51sq = getSinSqTerm(myMin.model.dm51Sq,E1,E2,l1,l2);
		sin54sq = getSinSqTerm(myMin.model.dm54Sq,E1,E2,l1,l2);
	}
	else{
		sin51sq = 0;
		sin54sq = 0;
	}

	if(myMin.model.dm61Sq != 0){
		sin61sq = getSinSqTerm(myMin.model.dm61Sq,E1,E2,l1,l2);
		sin64sq = getSinSqTerm(myMin.model.dm64Sq,E1,E2,l1,l2);
		sin65sq = getSinSqTerm(myMin.model.dm65Sq,E1,E2,l1,l2);
	}
	else{
		sin61sq = 0;
		sin64sq = 0;
		sin65sq = 0;
	}

	double osc_int = 4 * ((1 - pow(myMin.model.Ue[0],2) - pow(myMin.model.Ue[1],2) - pow(myMin.model.Ue[2],2))
			* (pow(myMin.model.Ue[0],2) * sin41sq + pow(myMin.model.Ue[1],2) * sin51sq + pow(myMin.model.Ue[2],2)*sin61sq)
			+ pow(myMin.model.Ue[0],2) * pow(myMin.model.Ue[1],2) * sin54sq
			+ pow(myMin.model.Ue[0],2) * pow(myMin.model.Ue[2],2) * sin64sq
			+ pow(myMin.model.Ue[1],2) * pow(myMin.model.Ue[2],2) * sin65sq);

	return osc_int;
}
void fcnXsec(int &npar, double *gin, double &fval, double *xval, int iflag){

	double *_gin = gin;	int &_npar = npar;	int _iflag = iflag;

	int nBins_lsnd = 5;
	int nBins_karmen = 6;
	double d_karmen = 16.2;	double l_karmen = 3.;
	double d_lsnd = 25.475;	double l_lsnd = 8.75;

	double k_lsnd  = xval[0];
	double k_karmen = xval[1];
	double k_correl = xval[2];

	double chisq;
	chisq = pow((k_lsnd - 1.)/myMin.xPack.lsnd_sys,2) + pow((k_karmen - 1.)/myMin.xPack.karmen_sys,2) + pow((k_correl - 1.)/myMin.xPack.correl_sys,2);

	double delE, lLo, lHi, ELo, EHi, EAvg, oscProb, xSec;
	// First, let's loop through karmen
	for(int i = 0; i < nBins_karmen; i++){
		if(i > 0)
			delE = (myMin.xPack.karmen_Enu[i] - myMin.xPack.karmen_Enu[i-1])/2.;
		else
			delE = (myMin.xPack.karmen_Enu[i+1] - myMin.xPack.karmen_Enu[i])/2.;
		lLo = d_karmen;
		lHi = d_karmen + l_karmen;
		ELo = myMin.xPack.karmen_Enu[i] - delE;
		EHi = myMin.xPack.karmen_Enu[i] + delE;
		EAvg = (EHi + ELo)/2.;
		myMin.xPack.ESigma = max(.08, .115/sqrt(EAvg - 17.3));
		oscProb = osc_int(ELo,EHi,lLo,lHi);
		xSec = (-2.5954e-4 * pow(EAvg,3) + 5.0028e-2 * pow(EAvg,2) - 1.5280 * EAvg + 12.876) * (1 - oscProb) * k_karmen * k_correl;
		chisq += pow((myMin.xPack.karmen[i] - xSec)/myMin.xPack.karmen_error[i],2);
	}

	// Now, do whatever we just did, but to LSND
	for(int i = 0; i < nBins_lsnd; i++){
		if(i > 0)
			delE = (myMin.xPack.lsnd_Enu[i] - myMin.xPack.lsnd_Enu[i-1])/2.;
		else
			delE = (myMin.xPack.lsnd_Enu[i+1] - myMin.xPack.lsnd_Enu[i])/2.;
		lLo = d_lsnd;
		lHi = d_lsnd + l_lsnd;
		ELo = myMin.xPack.lsnd_Enu[i] - delE;
		EHi = myMin.xPack.lsnd_Enu[i] + delE;
		EAvg = (EHi + ELo)/2.;
		myMin.xPack.ESigma = max(.08, .48/sqrt(EAvg - 17.3));
		oscProb = osc_int(ELo,EHi,lLo,lHi);
		xSec = (-2.5954e-4 * pow(EAvg,3) + 5.0028e-2 * pow(EAvg,2) - 1.5280 * EAvg + 12.876) * (1 - oscProb) * k_lsnd * k_correl;
		chisq += pow((myMin.xPack.lsnd[i] - xSec)/myMin.xPack.lsnd_error[i],2);
	}

	fval = chisq;

	return;
}
chisqStruct getChi2Xsec(neutrinoModel model, xsecPackage pack){

	chisqStruct result;
	result.zero();

	double chisq = 0;

	myMin.xPack = pack;
	myMin.model = model;

	double kkarmen, klsnd, kcorrel;

	double arglis[2];
    int ierflag;
    arglis[0] = -1.;
	gMinuit->SetFCN(fcnXsec);
    // Okay, let's get this minuit garbage started
    gMinuit->mnexcm("SET PRI",arglis,1,ierflag);
	gMinuit->mnexcm("SET NOWARNING",arglis,0,ierflag);
	arglis[0] = 1.;
	gMinuit->mnexcm("SET STR",arglis,1,ierflag);
    // Now, clear parameters and set 'em anew
    gMinuit->mnexcm("CLE",arglis,1,ierflag);
	gMinuit->mnparm(0,TString("k_lsnd"),1.,.1,0.,0.,ierflag);
	gMinuit->mnparm(1,TString("k_karmen"),1.,.1,0.,0.,ierflag);
	gMinuit->mnparm(2,TString("k_correl"),1.,.1,0.,0.,ierflag);

	// With these high dm2's, CHOOZ oscillations average out so we don't need to fit energy scale factor g
	arglis[0] = 5.;
	gMinuit->mnexcm("CALL", arglis,1,ierflag);
	double disttomin, errdef;   int npari, nparx, istat;
	arglis[0] = 100000.;
	arglis[1] = .1;
	gMinuit->mnexcm("MINI", arglis,1,ierflag);

	gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
	result.chi2 = chisq;

	return result;
}

// Bugey, finally working more or less
void fcnBugey(int &npar, double *gin, double &fval, double  *xval, int iflag){

	double *_gin = gin;	int &_npar = npar;	int _iflag = iflag;

  	double normReactorAno[] = {1.06237,1.06197,1.0627};
  	double nBins[] = {25, 25, 10};
  	double bigA, b, smallA[3];

  	bigA = xval[0];
  	smallA[0] = xval[1];
  	smallA[1] = xval[2];
  	smallA[2] = xval[3];
  	b = xval[4];

  	double prob, chisq;
  	double Ei = 1.;
  	double sinSq[dm2VecMaxDim];
  	ROOT::Math::Interpolator dif(dm2VecMaxDim,ROOT::Math::Interpolation::kCSPLINE);
	oscContribution oscCon = getOscContributionsNueDis(myMin.model);

  	chisq = 0.;

  	for(int j = 0; j < 3; j++){ // loop over baselines
    	for(int i = 0; i < nBins[j]; i++){  // loop over energy bins
      		for(int k = 0; k < dm2VecMaxDim; k++){
        		sinSq[k] = myMin.bPack.sinSqDeltaGrid[k][i][j];
      		}
      		// We're doing nue disappearance
			prob = 1.;
      		dif.SetData(dm2VecMaxDim,dm2Vec,sinSq);
      		for(int iCon = 0; iCon < 6; iCon ++){
				if(oscCon.dm2[iCon] != 0.){
					prob += oscCon.aEE[iCon] * dif.Eval(oscCon.dm2[iCon]);
				}
			}

      		// Now, calculate the chisq
      		double num = (bigA * smallA[j] + b * (myMin.bPack.energy[j][i] - Ei)) * prob - myMin.bPack.observed[j][i];
			if(ReactorAnomaly)
    			num = (bigA * smallA[j] + b * (myMin.bPack.energy[j][i] - Ei)) * prob * normReactorAno[j] - myMin.bPack.observed[j][i];
			chisq += pow(num/myMin.bPack.sigmaRatio[j][i],2);

    	}

    	chisq += pow((smallA[j] - 1.)/myMin.bPack.sigmaSmallA,2);
  	}

  	chisq += pow((bigA - 1.)/myMin.bPack.sigmaBigA,2) + pow(b/myMin.bPack.sigmaB,2);
  	fval = chisq;
}
chisqStruct getChi2Bugey(neutrinoModel model, bugeyPackage pack){

  	chisqStruct result;
  	result.zero();

  	double chisq = 0;

  	myMin.bPack = pack;
  	myMin.model = model;

  	double arglis[2];
  	int ierflag;
  	arglis[0] = -1.;
  	gMinuit->SetFCN(fcnBugey);
  	gMinuit->mnexcm("SET PRI",arglis,1,ierflag);
  	// Okay, let's get this minuit garbage started
  	arglis[0] = 1.;
  	gMinuit->mnexcm("SET STR",arglis,1,ierflag);
  	// Now, clear parameters and set 'em anew
  	gMinuit->mnexcm("CLE",arglis,1,ierflag);
  	gMinuit->mnparm(0,TString("bigA"),1.,pack.sigmaBigA,.5,1.5,ierflag);
  	gMinuit->mnparm(1,TString("smallA1"),1.,pack.sigmaSmallA,.5,1.5,ierflag);
  	gMinuit->mnparm(2,TString("smallA2"),1.,pack.sigmaSmallA,.5,1.5,ierflag);
  	gMinuit->mnparm(3,TString("smallA3"),1.,pack.sigmaSmallA,.5,1.5,ierflag);
  	gMinuit->mnparm(4,TString("b"),0.,pack.sigmaB,-.1,.1,ierflag);

  	// Now reset the function call and check if we even need to minimize. If it's a lost cause, don't waste the computation time.
  	arglis[0] = 5.;
  	gMinuit->mnexcm("CALL",arglis,1,ierflag);
  	double disttomin, errdef;   int npari, nparx, istat;

  	arglis[0] = 10000.;
  	arglis[1] = 100.;
  	gMinuit->mnexcm("MINI",arglis,2,ierflag);
  	gMinuit->mnstat(chisq,disttomin,errdef,npari,nparx,istat);
  	result.chi2 = chisq;

  	return result;
}
