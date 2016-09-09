/* ------------------------------------------//
Created by Davio Cianci and Georgia Karagiorgi
Jan 25th, 2016

------------------------------------------// */

// Declare JobOptions Variables
int noOfSteriles, CPConserving, scanType, gridPoints, jobID, nMCGen, rndInit, BugeyProcess, CCFRProcess, CDHSProcess, CHOOZProcess, KARMENProcess, LSNDProcess, NOMADProcess, MBProcess, MBProcessNubar, ATMOSPHERICProcess,
        NUMIProcess, MINOSProcess, MINOSNCProcess, GALLIUMProcess, ReactorAnomaly, XSECProcess, MBDISProcess, MBDISProcessNubar, nRuns;
double chi2Cut, stepSize, temperature, UMax, UMaxSq, trim;
bool usingUe = true;
bool usingUm = true;
std::string jobOptLoc = "inputs/";
std::string dataLoc = "data/";
bool debug = false;

// For Integral Evaluations,
double dm2Vec[dm2VecMaxDim];
void dm2VecInit(double dm2Min, double dm2Max){
    for(int i = 0; i < dm2VecMaxDim; i++){
        dm2Vec[i] = pow(10,TMath::Log10(dm2Min) + double(i) / (dm2VecMaxDim-1) * TMath::Log10(dm2Max/dm2Min));
    }
}

// Number of degrees of freedom
int ndf = 0;
void getNDF(){

	std::cout << "Total bins: " << ndf << std::endl;

    ndf -= noOfSteriles;
    if(BugeyProcess + CHOOZProcess + KARMENProcess + LSNDProcess +
            NOMADProcess + MBProcess + MBProcessNubar + NUMIProcess +
            GALLIUMProcess + XSECProcess == 0)
        usingUe = false;
	else	ndf -= noOfSteriles;
    if (CCFRProcess + CDHSProcess + KARMENProcess + LSNDProcess +
            NOMADProcess + MBProcess + MBProcessNubar + NUMIProcess +
            MINOSProcess+MINOSNCProcess+ATMOSPHERICProcess+
            MBDISProcess + MBDISProcessNubar == 0)
        usingUm = false;
	else	ndf -= noOfSteriles;
    if (BugeyProcess + CHOOZProcess + CCFRProcess + CDHSProcess + GALLIUMProcess + MINOSProcess +
            MINOSNCProcess + MBDISProcess + MBDISProcessNubar == 0)
        ndf += noOfSteriles;
    if(KARMENProcess + LSNDProcess + NOMADProcess + NUMIProcess + MBProcess + MBProcessNubar > 0 && CPConserving == 0){
        if(noOfSteriles == 2)   ndf -= 1;
        if(noOfSteriles == 3)   ndf -= 3;
    }
}

// Detector Initializations
boonePackage mbNuInit(){
    boonePackage pack;

    int nBins = 11;
    int nBins_mu = 8;
    pack.nFOscEvts = 17204;
    pack.full_fractCovMatrix.resize(nBins + nBins + nBins_mu, std::vector<float>(nBins + nBins + nBins_mu));
    pack.EnuQE = new float[nBins+1];
    pack.NueData = new int[nBins+nBins];
    pack.NumuData = new int[nBins_mu];
    pack.NueBgr = new float[nBins];
    pack.Numu = new float[nBins_mu];
    pack.FOsc_EnuQE = new float[pack.nFOscEvts];     // reconstructed neutrino energy
    pack.FOsc_EnuTrue = new float[pack.nFOscEvts];   // true energy of neutrino
    pack.FOsc_LnuTrue = new float[pack.nFOscEvts];   // distance from production and detection points
    pack.FOsc_weight = new float[pack.nFOscEvts];

    ifstream file;
    file.open(dataLoc+"miniboone_binboundaries_lowe.txt");
    for(int i = 0; i < nBins+1; i++)
        file >> pack.EnuQE[i];
    file.close();

    // Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
    file.open(dataLoc+"miniboone_nuedata_lowe.txt");
    for(int i = 0; i < nBins; i++)
        file >> pack.NueData[i];
    file.close();

    // Get measured Numu ccqe events per enuqe bin
    file.open(dataLoc+"miniboone_numudata.txt");
    for(int i = 0; i < nBins_mu; i++)
        file >> pack.NumuData[i];
    file.close();

    // Get predicted nue background events per enuqe bin
    file.open(dataLoc+"miniboone_nuebgr_lowe.txt");
    for(int i = 0; i < nBins; i++)
        file >> pack.NueBgr[i];
    file.close();

    // Get predicted numu ccqe events per enuqe bin
    file.open(dataLoc+"miniboone_numu.txt");
    for(int i = 0; i < nBins_mu; i++)
        file >> pack.Numu[i];
    file.close();

    // Get fractional cov matrix for full numu->nue oscillation events
    file.open(dataLoc+"neutrino_frac_error_matrix.txt");
    for(int i = 0; i < nBins + nBins + nBins_mu; i++)
        for(int j = 0; j < nBins + nBins + nBins_mu; j++)
            file >> pack.full_fractCovMatrix[i][j];
    file.close();

    // Get nFullOscEvts for full numu->nue osc events after nue cuts
    file.open(dataLoc+"miniboone_numode_fullosc_ntuple.txt");
    for(int i = 0; i < pack.nFOscEvts; i++){
        file >> pack.FOsc_EnuQE[i];     // reconstructed neutrino energy
        file >> pack.FOsc_EnuTrue[i];   // true energy of neutrino
        file >> pack.FOsc_LnuTrue[i];   // distance from production and detection points
        file >> pack.FOsc_weight[i];    // event weight
    }
    file.close();

    ndf += nBins + nBins_mu - 1;
	if(debug) std::cout << "MBnu initialized. Bins: " << nBins + nBins_mu - 1 << std::endl;
    return pack;
}
boonePlusPackage mbNuInitPlus(){
	boonePlusPackage pack;

	const short nBins = 11;
	const short nBins_mu = 8;
	const int nu_nfosc = 17204;
	pack.nFOscEvts = nu_nfosc;
	pack.full_fractCovMatrix.resize(nBins + nBins + nBins_mu, std::vector<float>(nBins + nBins + nBins_mu));
	pack.NueData = new int[nBins+nBins];
	pack.NumuData = new int[nBins_mu];
	pack.NueBgr = new float[nBins];
	pack.Numu = new float[nBins_mu];
	float *nu_EnuQE = new float[nBins + 1];
	float *nu_FOsc_EnuQE = new float[nu_nfosc];
	float *nu_FOsc_EnuTrue = new float[nu_nfosc];
	float *nu_FOsc_LnuTrue = new float[nu_nfosc];
	float *nu_FOsc_weight = new float[nu_nfosc];

	ifstream file;
	// Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
	file.open(dataLoc+"miniboone_nuedata_lowe.txt");
	for(int i = 0; i < nBins; i++)
		file >> pack.NueData[i];
	file.close();

	// Get measured Numu ccqe events per enuqe bin
	file.open(dataLoc+"miniboone_numudata.txt");
	for(int i = 0; i < nBins_mu; i++)
		file >> pack.NumuData[i];
	file.close();

	// Get predicted nue background events per enuqe bin
	file.open(dataLoc+"miniboone_nuebgr_lowe.txt");
	for(int i = 0; i < nBins; i++)
		file >> pack.NueBgr[i];
	file.close();

	// Get predicted numu ccqe events per enuqe bin
	file.open(dataLoc+"miniboone_numu.txt");
	for(int i = 0; i < nBins_mu; i++)
		file >> pack.Numu[i];
	file.close();

	// Get fractional cov matrix for full numu->nue oscillation events
	file.open(dataLoc+"neutrino_frac_error_matrix.txt");
	for(int i = 0; i < nBins + nBins + nBins_mu; i++)
		for(int j = 0; j < nBins + nBins + nBins_mu; j++)
			file >> pack.full_fractCovMatrix[i][j];
	file.close();

	file.open(dataLoc+"miniboone_binboundaries_lowe.txt");
	for(int i = 0; i < nBins+1; i++)
		file >> nu_EnuQE[i];
	file.close();

	// Get nFullOscEvts for full numu->nue osc events after nue cuts
	file.open(dataLoc+"miniboone_numode_fullosc_ntuple.txt");
	for(int iEvt = 0; iEvt < nu_nfosc; iEvt++){
        file >> nu_FOsc_EnuQE[iEvt];
        file >> nu_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nu_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nu_FOsc_weight[iEvt];    // event weight
	}
	file.close();

	float dmmax = 100.;
	float dmmin = 0.01;
	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	pack.lib_sinsq.resize(100, std::vector<float>(nBins));
	pack.lib_sin.resize(100, std::vector<float>(nBins));

	for(int mi = 0; mi < 100; mi++){
		std::cout << "mass numbah: " << mi << std::endl;

		float dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		for(int iB = 0; iB < nBins; iB++){
			pack.lib_sinsq[mi][iB] = 0;
			pack.lib_sin[mi][iB] = 0;
		}

		for(int iFOsc = 0; iFOsc < nu_nfosc; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	            if(nu_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nu_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){
					float ETru = nu_FOsc_EnuTrue[iFOsc];
					float LTru = nu_FOsc_LnuTrue[iFOsc];

					pack.lib_sinsq[mi][iB] += nu_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
					pack.lib_sin[mi][iB] += nu_FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
	            }
	        }
	    }
	}

	ndf += nBins + nBins_mu - 1;
	if(debug) std::cout << "MBnu initialized. Bins: " << nBins + nBins_mu - 1 << std::endl;
	return pack;
}
boonePackage mbNubarInit(){
    boonePackage pack;

    int nBins = 11;
    int nBins_mu = 8;
    pack.nFOscEvts = 86403;

    pack.full_fractCovMatrix.resize(nBins + nBins + nBins_mu, std::vector<float>(nBins + nBins + nBins_mu));
    pack.EnuQE = new float[nBins+1];
    pack.NueData = new int[nBins+nBins];
    pack.NumuData = new int[nBins_mu];
    pack.NueBgr = new float[nBins];
    pack.Numu = new float[nBins_mu];
    pack.FOsc_EnuQE = new float[pack.nFOscEvts];     // reconstructed neutrino energy
    pack.FOsc_EnuTrue = new float[pack.nFOscEvts];   // true energy of neutrino
    pack.FOsc_LnuTrue = new float[pack.nFOscEvts];   // distance from production and detection points
    pack.FOsc_weight = new float[pack.nFOscEvts];

    ifstream file;
    file.open(dataLoc+"miniboone_binboundaries_lowe.txt");
    for(int i = 0; i < nBins+1; i++)
        file >> pack.EnuQE[i];
    file.close();

    // Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
    file.open(dataLoc+"miniboone_nuebardata_lowe.txt");
    for(int i = 0; i < nBins; i++)
        file >> pack.NueData[i];
    file.close();

    // Get measured Numu ccqe events per enuqe bin
    file.open(dataLoc+"miniboone_numubardata.txt");
    for(int i = 0; i < nBins; i++)
        file >> pack.NumuData[i];
    file.close();

    // Get predicted nue background events per enuqe bin
    file.open(dataLoc+"miniboone_nuebarbgr_lowe.txt");
    for(int i = 0; i < nBins; i++)
        file >> pack.NueBgr[i];
    file.close();

    // Get predicted numu ccqe events per enuqe bin
    file.open(dataLoc+"miniboone_numubar.txt");
    for(int i = 0; i < nBins_mu; i++)
        file >> pack.Numu[i];
    file.close();

    // Get fractional cov matrix for full numu->nue oscillation events
    file.open(dataLoc+"miniboone_full_fractcovmatrix_nubar_lowe.txt");
    for(int i = 0; i < nBins + nBins + nBins_mu; i++)
        for(int j = 0; j < nBins + nBins + nBins_mu; j++)
            file >> pack.full_fractCovMatrix[i][j];
    file.close();

    // Get nFullOscEvts for full numu->nue osc events after nue cuts
    file.open(dataLoc+"miniboone_numubarnuebarfullosc_ntuple.txt");
    for(int i = 0; i < pack.nFOscEvts; i++){
        file >> pack.FOsc_EnuQE[i];     // reconstructed neutrino energy
        file >> pack.FOsc_EnuTrue[i];   // true energy of neutrino
        file >> pack.FOsc_LnuTrue[i];   // distance from production and detection points
        file >> pack.FOsc_weight[i];    // event weight
    }
    file.close();

    ndf += nBins + nBins_mu - 1;
	if(debug) std::cout << "MBnubar initialized. Bins: " << nBins + nBins_mu - 1 << std::endl;
    return pack;
}
boonePlusPackage mbNubarInitPlus(){
	boonePlusPackage pack;

	const short nBins = 11;
	const short nBins_mu = 8;
	const int nubar_nfosc = 86403;
	pack.nFOscEvts = nubar_nfosc;
	pack.full_fractCovMatrix.resize(nBins + nBins + nBins_mu, std::vector<float>(nBins + nBins + nBins_mu));
	pack.NueData = new int[nBins+nBins];
	pack.NumuData = new int[nBins_mu];
	pack.NueBgr = new float[nBins];
	pack.Numu = new float[nBins_mu];
	float *nubar_EnuQE = new float[nBins + 1];
	float *nubar_FOsc_EnuQE = new float[nubar_nfosc];
	float *nubar_FOsc_EnuTrue = new float[nubar_nfosc];
	float *nubar_FOsc_LnuTrue = new float[nubar_nfosc];
	float *nubar_FOsc_weight = new float[nubar_nfosc];

	ifstream file;
	// Get measured Nue events per reconstructed electron neutrino energy bin (enuqe)
	file.open(dataLoc+"miniboone_nuebardata_lowe.txt");
	for(int i = 0; i < nBins; i++)
		file >> pack.NueData[i];
	file.close();

	// Get measured Numu ccqe events per enuqe bin
	file.open(dataLoc+"miniboone_numubardata.txt");
	for(int i = 0; i < nBins_mu; i++)
		file >> pack.NumuData[i];
	file.close();

	// Get predicted nue background events per enuqe bin
	file.open(dataLoc+"miniboone_nuebarbgr_lowe.txt");
	for(int i = 0; i < nBins; i++)
		file >> pack.NueBgr[i];
	file.close();

	// Get predicted numu ccqe events per enuqe bin
	file.open(dataLoc+"miniboone_numubar.txt");
	for(int i = 0; i < nBins_mu; i++)
		file >> pack.Numu[i];
	file.close();

	// Get fractional cov matrix for full numu->nue oscillation events
	file.open(dataLoc+"miniboone_full_fractcovmatrix_nubar_lowe.txt");
	for(int i = 0; i < nBins + nBins + nBins_mu; i++)
		for(int j = 0; j < nBins + nBins + nBins_mu; j++)
			file >> pack.full_fractCovMatrix[i][j];
	file.close();

	file.open(dataLoc+"miniboone_binboundaries_lowe.txt");
	for(int i = 0; i < nBins+1; i++)
		file >> nubar_EnuQE[i];
	file.close();

	// Get nFullOscEvts for full numu->nue osc events after nue cuts
	file.open(dataLoc+"miniboone_numubarnuebarfullosc_ntuple.txt");
	for(int iEvt = 0; iEvt < nubar_nfosc; iEvt++){
        file >> nubar_FOsc_EnuQE[iEvt];
        file >> nubar_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nubar_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nubar_FOsc_weight[iEvt];    // event weight
	}
	file.close();

	float dmmax = 100.;
	float dmmin = 0.01;
	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	pack.lib_sinsq.resize(100, std::vector<float>(nBins));
	pack.lib_sin.resize(100, std::vector<float>(nBins));

	for(int mi = 0; mi < 100; mi++){
		std::cout << "mass numbah: " << mi << std::endl;

		float dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		for(int iB = 0; iB < nBins; iB++){
			pack.lib_sinsq[mi][iB] = 0;
			pack.lib_sin[mi][iB] = 0;
		}

		for(int iFOsc = 0; iFOsc < nubar_nfosc; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	            if(nubar_FOsc_EnuQE[iFOsc] > nubar_EnuQE[iB] && nubar_FOsc_EnuQE[iFOsc] < nubar_EnuQE[iB+1]){
					float ETru = nubar_FOsc_EnuTrue[iFOsc];
					float LTru = nubar_FOsc_LnuTrue[iFOsc];

					pack.lib_sinsq[mi][iB] += nubar_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru*.01/ETru),2);
					pack.lib_sin[mi][iB] += nubar_FOsc_weight[iFOsc]*sin(1.267*2*dm2*LTru*.01/ETru);
	            }
	        }
	    }
	}

	ndf += nBins + nBins_mu - 1;
	if(debug) std::cout << "MBnubar initialized. Bins: " << nBins + nBins_mu - 1 << std::endl;
	return pack;
}
atmPackage atmInit(){
    atmPackage pack;

    double dmuMin = 0.;
    double dmuMax = .25;
    const int dmuVecMaxDim = 101;

    double temp0[] = {0.00000,0.03913,0.10552,0.18864,0.28771,0.40839,0.54925,0.69901,0.86038,1.05169,1.27231,1.53027,1.77613,2.02014,2.27515,2.54545,2.83366,3.13347,3.44052,3.76300,
                4.09784,4.46017,4.84822,5.26111,5.69840,6.16382,6.65471,7.16706,7.70432,8.24987,8.80458,9.38657,9.99066,10.61653,11.26452,11.93572,12.62999,13.32833,14.03633,14.77024,
                15.52484,16.29127,17.07689,17.87941,18.69526,19.53419,20.39532,21.27474,22.17658,23.08987,24.01778,24.969511,25.93752,26.92796,27.94275,28.98161,30.04210,31.12363,32.22914,33.35148,
                34.49207,35.65003,36.82462,38.01945,39.23037,40.46261,41.71198,42.97714,44.26171,45.56482,46.88852,48.22036,49.56687,50.93319,52.30869,53.58208,54.86813,56.16923,57.48430,58.81740,
                60.16430,61.52579,62.90204,64.28466,65.67754,67.07690,68.48444,69.90443,71.33371,72.77951,74.23951,75.70818,77.18138,78.66365,80.15671,81.66008,83.17423,84.69586,86.22622,87.77160,89.32895};
    pack.dmuVec = new double[dmuVecMaxDim];
    pack.dchi2Vec = new double[dmuVecMaxDim];

    for(int iA = 0; iA < dmuVecMaxDim; iA++){
        pack.dmuVec[iA] = dmuMin + (dmuMax - dmuMin) * (iA)/(dmuVecMaxDim-1);
        pack.dchi2Vec[iA] = temp0[iA];
    }

    ndf += 1;
	if(debug) std::cout << "Atmospheric initialized. Bins: " << 1 << std::endl;
    return pack;
}
numiPackage numiInit(){
    numiPackage pack;

    const int nBins = 10;
    const int nFOscEvts = 3323;

    pack.FOsc_EnuQE = new double[nFOscEvts];     // reconstructed neutrino energy
    pack.FOsc_EnuTrue = new double[nFOscEvts];   // true energy of neutrino
    pack.FOsc_LnuTrue = new double[nFOscEvts];   // distance from production and detection points
    pack.FOsc_weight = new double[nFOscEvts];
    pack.FOsc_fracError = new double[nBins];
	pack.EnuQE = new double[nBins+1];
	pack.NueData = new int[nBins];
	pack.NueBgr = new double[nBins];
	pack.NueBgr_error = new double[nBins];

    double temp0[] = {200.,300.,475.,675.,900.,1100.,1300.,1500.,1700.,2000.,3000.};
    int temp1[] = {59,142,151,146,83,68,57,39,19,16};
    double temp2[] = {41.7574,117.537,123.673,118.025,82.0508,64.8166,41.1826,31.3267,22.0301,17.6672};
    double temp3[] = {12.112,25.5363,27.4975,24.8745,18.1301,15.11,10.5605,9.03036,6.88422,5.70684};

    for(int i = 0; i < nBins; i++){
    	pack.EnuQE[i] = temp0[i];
    	pack.NueData[i] = temp1[i];
    	pack.NueBgr[i] = temp2[i];
    	pack.NueBgr_error[i] = temp3[i];
    }
    pack.EnuQE[nBins] = temp0[nBins];

    ifstream file;
    // Get nFullOscEvts from another file!
    file.open(dataLoc+"numi_fullosc.out");
    for(int i = 0; i < nFOscEvts; i++){
        file >> pack.FOsc_weight[i];     // reconstructed neutrino energy
        file >> pack.FOsc_EnuTrue[i];   // true energy of neutrino
        file >> pack.FOsc_EnuQE[i];   // distance from production and detection points
        file >> pack.FOsc_LnuTrue[i];    // event weight
    }
    file.close();

    for(int i = 0; i < nBins; i++){
        pack.FOsc_fracError[i] = pack.NueBgr_error[i] / pack.NueBgr[i];
    }

    ndf += nBins;
	if(debug) std::cout << "Numi initialized. Bins: " << nBins << std::endl;
    return pack;
}
sinSqPackage lsndInit(){
    sinSqPackage pack;

    const int nBins = 5;
    double EMin = 20.; double EMax = 60.;

	pack.dm2Vec.resize(601);
    pack.sinSqDeltaGrid.resize(dm2VecMaxDim, std::vector<double>(nBins));
    pack.sinSqDeltaGrid2.resize(dm2VecMaxDim, std::vector<double>(nBins));

    double temp0[] = {10., 11.7, 17.6, 7.8, 3.7};
    double temp1[] = {5.1, 5.2, 3.7, 2.2, 0.6};
    pack.observed = new double[nBins];
    pack.bkg = new double[nBins];
    for(int iL = 0; iL < nBins; iL++){
    	pack.observed[iL] = temp0[iL];
    	pack.bkg[iL] = temp1[iL];
    }

    double energyMin[nBins], energyMax[nBins];
    for(int iL = 0; iL < nBins; iL++){
        energyMin[iL] = EMin + (double(iL)/double(nBins)) * (EMax - EMin);
        energyMax[iL] = EMin + (double(iL+1)/double(nBins)) * (EMax - EMin);
    }

    // Now, let's find the LSND integrals!
    double sum = 0.;
    double fullConversionLSNDTotal = 33300.;
    double rGammaCutEfficiency = .39;

    integralFuncsLSND funcLSND;
	/*
    for(int iL = 0; iL < nBins; iL++){
        ROOT::Math::Functor LSNDNorm(&funcLSND, &integralFuncsLSND::normFunc,3);
        double a[3] = {energyMin[iL],20,30. - 3.5};
        double b[3] = {energyMax[iL],60,30. + 3.5};

        ROOT::Math::IntegratorMultiDim ig;
        ig.SetFunction(LSNDNorm);
        sum += ig.Integral(a,b);
    }
    pack.norm = fullConversionLSNDTotal * rGammaCutEfficiency / sum;
	*/
	pack.norm = 57401.252; // Taken from old code since multiple integrals in root are tricky and time consuming

    /*
    for(int k = 0; k < dm2VecMaxDim; k++){
        for(int iL = 0; iL < nBins; iL++){

            funcLSND._dm2 = dm2Vec[k];

            ROOT::Math::Functor LSNDSinSq(&funcLSND, &integralFuncsLSND::sinSqFunction,3);
            ROOT::Math::Functor LSNDSinSq2(&funcLSND, &integralFuncsLSND::sinSqFunctionCPV,3);

            double a[] = {energyMin[iL],20,30. - 3.5};
            double b[] = {energyMax[iL],60,30. + 3.5};
            ROOT::Math::IntegratorMultiDim ig;
            ig.SetFunction(LSNDSinSq);
            pack.sinSqDeltaGrid[k][iL] = ig.Integral(a,b);// / (energyMax[iL] - energyMin[iL]);
            ig.SetFunction(LSNDSinSq2);
            pack.sinSqDeltaGrid2[k][iL] = ig.Integral(a,b);// / (energyMax[iL] - energyMin[iL]);

			if(k == 0) {
				//std::cout << energyMin[iL] << " " << energyMax[iL] << std::endl;
				std::cout << pack.sinSqDeltaGrid[k][iL] << " " << pack.sinSqDeltaGrid2[k][iL] << std::endl;
			}
		}
    }
	*/
	ifstream file;
	file.open(dataLoc+"lsndsinsq.txt");
	for(int k = 0; k < 601; k++){
        for(int iL = 0; iL < nBins; iL++){
			file >> pack.dm2Vec[k];
			file >> pack.sinSqDeltaGrid[k][iL];
			file >> pack.sinSqDeltaGrid2[k][iL];
		}
	}
	file.close();

    ndf += nBins;
	if(debug) std::cout << "LSND initialized. Bins: " << nBins << std::endl;
    return pack;
}
sinSqPackage karmenInit(){
    sinSqPackage pack;

    const int nBins = 9;

    pack.sinSqDeltaGrid.resize(dm2VecMaxDim, std::vector<double>(nBins));
    pack.sinSqDeltaGrid2.resize(dm2VecMaxDim, std::vector<double>(nBins));

    double temp0[] = {3.,4.,1.,3.,3.,1.,0.,0.,0.};
    double temp1[] = {3.4,3.7,3.4,2.3,1.1,0.8,0.6,0.4,0.1};
    pack.observed = new double[nBins];
    pack.bkg = new double[nBins];
    for(int iL = 0; iL < nBins; iL++){
    	pack.observed[iL] = temp0[iL];
    	pack.bkg[iL] = temp1[iL];
    }

    double EMin = 16.; double EMax = 52.;
    double energyMin[nBins], energyMax[nBins];
    for(int iK = 0; iK < nBins; iK++){
        energyMin[iK] = EMin + (double(iK)/double(nBins)) * (EMax - EMin);
        energyMax[iK] = EMin + (double(iK + 1)/double(nBins)) * (EMax - EMin);
    }

    // Now, let's find the KARMEN integrals!
    double fullMixHighDm2KarmenTotal = 2913.;
    double sum = 0.;
    integralFuncsKarmen integralFuncs;

    for(int iK = 0; iK < nBins; iK++){
        EMin = energyMin[iK];
        EMax = energyMax[iK];

        ROOT::Math::Functor1D wf(&integralFuncs,&integralFuncsKarmen::normFunc);
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
        ig.SetAbsTolerance(1*pow(10,-4));
        ig.SetRelTolerance(1*pow(10,-4));
        ig.SetFunction(wf);
        sum += ig.Integral(EMin,EMax);
    }
    pack.norm = fullMixHighDm2KarmenTotal / sum;

    for(int k = 0; k < dm2VecMaxDim; k++){
        for(int iK = 0; iK < nBins; iK++){
            integralFuncs._dm2 = dm2Vec[k];
            EMin = energyMin[iK];
            EMax = energyMax[iK];

            ROOT::Math::Functor1D wf1(&integralFuncs, &integralFuncsKarmen::sinSqFunction);
            ROOT::Math::Functor1D wf2(&integralFuncs, &integralFuncsKarmen::sinSqFunctionCPV);
            ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
            ig.SetAbsTolerance(1*pow(10,-4));
            ig.SetRelTolerance(1*pow(10,-4));
            ig.SetFunction(wf1);
            pack.sinSqDeltaGrid[k][iK] = ig.Integral(EMin, EMax);
            ig.SetFunction(wf2);
            pack.sinSqDeltaGrid2[k][iK] = ig.Integral(EMin, EMax);
        }
    }
    ndf += nBins;
	if(debug) std::cout << "Karmen initialized. Bins: " << nBins << std::endl;
    return pack;
}
galPackage galInit(){
    galPackage pack;

    // Gallium's nice and easy up here.
    double temp0[] = {1.,0.81,0.95,0.79};
    double temp1[] = {0.1,0.1,0.12,0.1};
    double temp2[] = {0.811,0.813};
    double temp3[] = {0.902,0.098};
    double temp4[] = {70.1,70.3};
    double temp5[] = {0.747,0.742,0.427,0.432};
    double temp6[] = {0.8163,0.0849,0.0895,0.0093};
    double temp7[] = {60.8,61.5,26.7,27.1};

    pack.obsRatioGal = new double[4];
    pack.errorGal = new double[4];
    pack.arLinesE = new double[2];
    pack.arLinesBr = new double[2];
    pack.arLinesXSec = new double[2];
    pack.crLinesE = new double[4];
    pack.crLinesBr = new double[4];
    pack.crLinesXSec = new double[4];

    for(int i = 0; i < 4; i++){
    	pack.obsRatioGal[i] = temp0[i];
    	pack.errorGal[i] = temp1[i];
    	pack.crLinesE[i] = temp5[i];
    	pack.crLinesBr[i] = temp6[i];
    	pack.crLinesXSec[i] = temp7[i];
    }
    for(int i = 0; i < 2; i++){
    	pack.arLinesE[i] = temp2[i];
    	pack.arLinesBr[i] = temp3[i];
    	pack.arLinesXSec[i] = temp4[i];
    }

	// Here's a little experiment.
	// so, L will go from 0 to the corner. In this way, we'll magically reduce two parameters to only one!
	// we need three L's and their corresponding volume integrals
	double nLength, nInt;
	nLength = 2000; nInt = 600;

	pack.volInt.resize(3, std::vector<double>(nLength));
	for(int i = 0; i < nLength; i++){
		pack.volInt[0][i] = 0.;	pack.volInt[1][i] = 0.;	pack.volInt[2][i] = 0.;
	}
	double radiusGallex = 1.9;  double heightGallex = 5.0;  double sourceHeightGallex[2] = {2.7,2.38};
	double radiusSage = 0.7;    double heightSage = 1.47;   double sourceHeightSage = 0.72;

	// Let's just start this whole mess over and make sure it works from the start.
	double Radius[3] = {radiusGallex, radiusGallex, radiusSage};
	double Height[3] = {heightGallex, heightGallex, heightSage};
	double SourceHeight[3] = {sourceHeightGallex[0], sourceHeightGallex[1], sourceHeightSage};

	for(int iex = 0; iex < 3; iex++){
		double maxht = max(Height[iex]-SourceHeight[iex],SourceHeight[iex]);
		double maxlen = sqrt(pow(Radius[iex],2) + pow(maxht,2));

		double ht, dht, rad, dr;
		for(int iHt = 0; iHt < nInt; iHt++){
			ht = (Height[iex] / nInt) * (iHt + .5); dht = Height[iex] / nInt;

			for(int iRd = 0; iRd < nInt; iRd++){
				rad = (Radius[iex] / nInt) * iRd; dr = Radius[iex] / nInt;

				int _ind = floor(sqrt(pow(rad,2) + pow(ht-SourceHeight[iex],2))/(maxlen/float(nLength)));
				pack.volInt[iex][_ind] += dr * dht * 2*TMath::Pi()*rad;
			}
		}
	}

    ndf += 4;
	if(debug) std::cout << "Gal initialized. Bins: " << 4 << std::endl;
    return pack;
}
minosPackage minosInit(){
    minosPackage pack;

    const int nBins = 12;
    const int nBins_ws = 13;

    pack.EnuQE = new double[nBins+1];
	pack.NumubarData = new double[nBins];
	pack.NumubarBkg = new double[nBins];
	pack.fracError = new double[nBins];
	pack.dataErr = new double[nBins];
	pack.EnuQE_ws = new double[nBins_ws+1];
	pack.NumubarData_ws = new double[nBins_ws];
	pack.NumubarBkg_ws = new double[nBins_ws];
	pack.fracError_ws = new double[nBins_ws];
	pack.dataErr_ws = new double[nBins_ws];

    // Minos looks easy, but DON'T BE FOOLED.
    // Right sign
    double temp0[] = {0, 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 15000., 20000.};
    double temp1[] = {0.14, 6.0, 13.1, 41.0, 29.0, 19.1, 7.9, 11.1, 11.1, 4.0, 4.8, 4.4};
    double temp2[] = {1.2, 20.9, 50.5, 56.8, 34.5, 17.8, 11.8, 9.1, 7.8, 7.0, 5.0, 2.7};
    double temp3[] = {0.020, 0.020, 0.06, 0.05, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.020, 0.020};
    double temp4[] = {1.04, 3.21, 4.30, 7.02, 6.10, 5.12, 3.54, 4.03, 4.03, 2.78, 1.31, 1.25};
    for(int i = 0; i < nBins; i++){
    	pack.EnuQE[i] = temp0[i];
    	pack.NumubarData[i] = temp1[i];
    	pack.NumubarBkg[i] = temp2[i];
    	pack.fracError[i] = temp3[i];
    	pack.dataErr[i] = temp4[i];
    }
	pack.EnuQE[nBins] = temp0[nBins];

    // Wrong sign
    double temp5[] = {0, 2000., 4000., 6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000., 25000., 30000., 50000.};
    double temp6[] = {0.0, 6.99, 22.00, 16.01, 14.95, 14.94, 13.94, 9.99, 7.94, 3.89, 1.94, 1.93, 0.73};
    double temp7[] = {3.4, 12.72, 16.31, 18.15, 18.99, 16.44, 13.94, 11.14, 8.84, 6.93, 4.38, 2.08, 0.53};
    double temp8[] = {.17, .13, .11, .1, .09, .09, .09, .11, .1, .11, .09, .16, .3};
    double temp9[] = {1.93, 3.29, 5.36, 4.67, 4.54, 4.52, 4.42, 3.79, 3.52, 2.67, 1.25, 1.2, .52};
    for(int i = 0; i < nBins_ws; i++){
    	pack.EnuQE_ws[i] = temp5[i];
    	pack.NumubarData_ws[i] = temp6[i];
    	pack.NumubarBkg_ws[i] = temp7[i];
    	pack.fracError_ws[i] = temp8[i];
    	pack.dataErr_ws[i] = temp9[i];
    }
	pack.EnuQE_ws[nBins_ws] = temp5[nBins_ws];

    ndf += nBins;
	if(debug) std::cout << "Minos initialized. Bins: " << nBins << std::endl;
    return pack;
}
minosncPackage minosncInit(){
    minosncPackage pack;

    pack.theta24_bestfit = 0.;
    pack.sd_theta24 = 0.087;

    ndf += 1;
    return pack;
}
booneDisPackage mbNuDisInit(){
    booneDisPackage pack;

	const short nBins = 16;
	const int nfosc = 1267007;
	pack.nFOscEvts = nfosc;

	pack.full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));
	pack.EnuQE = new float[nBins + 1];
	pack.NumuData = new float[nBins];
	pack.FOsc_EnuQE = new float[nfosc];
	pack.FOsc_EnuTrue = new float[nfosc];
	pack.FOsc_LnuTrue = new float[nfosc];
	pack.FOsc_weight = new float[nfosc];

	ifstream file;
	file.open(dataLoc+"miniboone_binboundaries_disap.txt");
	for(short i = 0; i < nBins+1; i++)
	    file >> pack.EnuQE[i];
	file.close();

	file.open(dataLoc+"miniboone_numudata_disap.txt");
	for(short i = 0; i < nBins; i++)
	    file >> pack.NumuData[i];
	file.close();

	file.open(dataLoc+"miniboone_frac_shape_matrix_numu_disap.txt");
	for(short i = 0; i < nBins; i++)
	    for(short j = 0; j < nBins; j++)
	        file >> pack.full_fractCovMatrix[i][j];
	file.close();

	file.open(dataLoc+"numudisap_ntuple.txt");
	int dummy;
	for(int iEvt = 0; iEvt < nfosc; iEvt++){
		file >> dummy;
	    file >> pack.FOsc_EnuQE[iEvt];
	    file >> pack.FOsc_EnuTrue[iEvt];   // true energy of neutrino
	    file >> pack.FOsc_LnuTrue[iEvt];   // distance from production and detection points
	    file >> pack.FOsc_weight[iEvt];    // event weight
	}
	file.close();

	ndf += nBins;
	if(debug) std::cout << "MBnu Dis initialized. Bins: " << nBins << std::endl;
	return pack;
}
booneDisPlusPackage mbNuDisInitPlus(){
	booneDisPlusPackage pack;

	const short nBins = 16;
	const int nu_nfosc = 1267007;

	pack.full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));
	pack.NumuData = new float[nBins];
    float *nu_EnuQE = new float[nBins + 1];
	float *nu_FOsc_EnuQE = new float[nu_nfosc];
	float *nu_FOsc_EnuTrue = new float[nu_nfosc];
	float *nu_FOsc_LnuTrue = new float[nu_nfosc];
	float *nu_FOsc_weight = new float[nu_nfosc];

	ifstream file;
	file.open(dataLoc+"miniboone_numudata_disap.txt");
	for(short i = 0; i < nBins; i++)
		file >> pack.NumuData[i];
	file.close();
	file.open(dataLoc+"miniboone_frac_shape_matrix_numu_disap.txt");
	for(short i = 0; i < nBins; i++)
		for(short j = 0; j < nBins; j++)
			file >> pack.full_fractCovMatrix[i][j];
	file.close();

    file.open(dataLoc+"miniboone_binboundaries_disap.txt");
    for(short i = 0; i < nBins+1; i++)
        file >> nu_EnuQE[i];
    file.close();

 	file.open(dataLoc+"numudisap_ntuple.txt");
	int dummy;
    for(int iEvt = 0; iEvt < nu_nfosc; iEvt++){
 		file >> dummy;
        file >> nu_FOsc_EnuQE[iEvt];
        file >> nu_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nu_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nu_FOsc_weight[iEvt];    // event weight
	}
	file.close();

	float dmmax = 100.;
	float dmmin = 0.01;
	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	pack.libdis_sinsq.resize(100, std::vector<float>(nBins));
	pack.libdis_noosc.resize(nBins);

	for(int mi = 0; mi < 100; mi++){
		std::cout << "mass numbah: " << mi << std::endl;

		float dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		for(int iB = 0; iB < nBins; iB++){
			pack.libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				pack.libdis_noosc[iB] = 0;
		}

		for(int iFOsc = 0; iFOsc < nu_nfosc; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	            if(nu_FOsc_EnuQE[iFOsc] > nu_EnuQE[iB] && nu_FOsc_EnuQE[iFOsc] < nu_EnuQE[iB+1]){

					float ETru = nu_FOsc_EnuTrue[iFOsc];
					float LTru = nu_FOsc_LnuTrue[iFOsc];

					pack.libdis_sinsq[mi][iB] += nu_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru/ETru),2);


					if(mi == 0)
						pack.libdis_noosc[iB] += nu_FOsc_weight[iFOsc];
	            }
	        }
	    }
	}
	return pack;
}
booneDisPackage mbNubarDisInit(){
    booneDisPackage pack;

    const short nBins = 16;
	const int nfosc = 686529;
	pack.nFOscEvts = nfosc;

    pack.full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));
    pack.EnuQE = new float[nBins + 1];
    pack.NumuData = new float[nBins];
	pack.FOsc_EnuQE = new float[nfosc];
	pack.FOsc_EnuTrue = new float[nfosc];
	pack.FOsc_LnuTrue = new float[nfosc];
	pack.FOsc_weight = new float[nfosc];

    ifstream file;
    file.open(dataLoc+"miniboone_binboundaries_disap.txt");
    for(short i = 0; i < nBins+1; i++)
        file >> pack.EnuQE[i];
    file.close();

    file.open(dataLoc+"miniboone_numubardata_disap.txt");
    for(short i = 0; i < nBins; i++)
        file >> pack.NumuData[i];
    file.close();

    file.open(dataLoc+"miniboone_frac_shape_matrix_numubar_disap.txt");
    for(short i = 0; i < nBins; i++)
        for(short j = 0; j < nBins; j++)
            file >> pack.full_fractCovMatrix[i][j];
    file.close();

	file.open(dataLoc+"numubardisap_ntuple.txt");
 	int dummy;
    for(int iEvt = 0; iEvt < nfosc; iEvt++){
 		file >> dummy;
        file >> pack.FOsc_EnuQE[iEvt];
        file >> pack.FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> pack.FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> pack.FOsc_weight[iEvt];    // event weight
	}
	file.close();

    ndf += nBins;
	if(debug) std::cout << "MBnubar Dis initialized. Bins: " << nBins << std::endl;
    return pack;
}
booneDisPlusPackage mbNubarDisInitPlus(){
	booneDisPlusPackage pack;

	const short nBins = 16;
	const int nfosc = 686529;

	pack.full_fractCovMatrix.resize(nBins, std::vector<float>(nBins));
	pack.NumuData = new float[nBins];
    float *nubar_EnuQE = new float[nBins + 1];
	float *nubar_FOsc_EnuQE = new float[nfosc];
	float *nubar_FOsc_EnuTrue = new float[nfosc];
	float *nubar_FOsc_LnuTrue = new float[nfosc];
	float *nubar_FOsc_weight = new float[nfosc];

	ifstream file;
	file.open(dataLoc+"miniboone_numubardata_disap.txt");
	for(short i = 0; i < nBins; i++)
		file >> pack.NumuData[i];
	file.close();
	file.open(dataLoc+"miniboone_frac_shape_matrix_numubar_disap.txt");
	for(short i = 0; i < nBins; i++)
		for(short j = 0; j < nBins; j++)
			file >> pack.full_fractCovMatrix[i][j];
	file.close();

    file.open(dataLoc+"miniboone_binboundaries_disap.txt");
    for(short i = 0; i < nBins+1; i++)
        file >> nubar_EnuQE[i];
    file.close();

 	file.open(dataLoc+"numubardisap_ntuple.txt");
	int dummy;
    for(int iEvt = 0; iEvt < nfosc; iEvt++){
 		file >> dummy;
        file >> nubar_FOsc_EnuQE[iEvt];
        file >> nubar_FOsc_EnuTrue[iEvt];   // true energy of neutrino
        file >> nubar_FOsc_LnuTrue[iEvt];   // distance from production and detection points
        file >> nubar_FOsc_weight[iEvt];    // event weight
	}
	file.close();

	float dmmax = 100.;
	float dmmin = 0.01;
	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	pack.libdis_sinsq.resize(100, std::vector<float>(nBins));
	pack.libdis_noosc.resize(nBins);

	for(int mi = 0; mi < 100; mi++){
		std::cout << "mass numbah: " << mi << std::endl;

		float dm2 = pow(10,((mi+1.)/100.*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		for(int iB = 0; iB < nBins; iB++){
			pack.libdis_sinsq[mi][iB] = 0;
			if(mi == 0)
				pack.libdis_noosc[iB] = 0;
		}

		for(int iFOsc = 0; iFOsc < nfosc; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred

	            if(nubar_FOsc_EnuQE[iFOsc] > nubar_EnuQE[iB] && nubar_FOsc_EnuQE[iFOsc] < nubar_EnuQE[iB+1]){

					float ETru = nubar_FOsc_EnuTrue[iFOsc];
					float LTru = nubar_FOsc_LnuTrue[iFOsc];

					pack.libdis_sinsq[mi][iB] += nubar_FOsc_weight[iFOsc]*pow(sin(1.267*dm2*LTru/ETru),2);
					if(mi == 0)
						pack.libdis_noosc[iB] += nubar_FOsc_weight[iFOsc];
	            }
	        }
	    }
	}
	return pack;
}
nomadPackage nomadInit(){
	nomadPackage pack;

    // Since everything's hard coded in, it's nice and easy for you, but a nightmare for me.
    const int maxEnergyBins = 10; const int maxRadialBins = 3;

    integralFuncsNomad integralFuncs;

    pack.sinSqDeltaGrid.resize(dm2VecMaxDim, std::vector<double>(30));
    pack.sinSqDeltaGrid2.resize(dm2VecMaxDim, std::vector<double>(30));
    pack.sigmaRemu.resize(30, std::vector<double>(30));

    double temp0[] = {0.00547,0.00446,0.00500,0.00523,0.00530,0.00592,0.00555,0.01111,0.01887,0.01051,0.00401,0.00454,0.00459,0.00545,0.00904,0.01519,0.01932,0.02306,0.02040,0.01002,0.00511,0.00488,0.00759,
                    0.01364,0.02365,0.03866,0.03610,0.03885,0.02489,0.01131};
	double temp1[] = {0.00516,0.00467,0.00493,0.00502,0.00489,0.00555,0.00660,0.01007,0.02049,0.01172,0.00489,0.00382,0.00492,0.00561,0.00846,0.01359,0.01899,0.02468,0.02243,0.01126,0.00463,0.00521,0.00724,0.01203,
                    0.02187,0.03343,0.03990,0.03893,0.02270,0.01092};
    pack.observed = new double[30];
	pack.bkg = new double[30];
    for(int i = 0; i < maxEnergyBins * maxRadialBins; i++){
    	pack.observed[i] = temp0[i];
    	pack.bkg[i] = temp1[i];
    }

    double sigmaSystRemu[30][30] = {
        {0.0150,0.0043,0.0045,0.0046,0.0045,0.0051,0.0060,0.0092,0.0187,0.0107,0.0142,0.0035,0.0045,0.0051,0.0077,0.0124,0.0173,0.0225,0.0204,0.0103,0.0134,0.0047,0.0066,0.0110,0.0199,0.0305,0.0364,0.0355,0.0207,0.0100},
        {0.0043,0.0091,0.0041,0.0041,0.0040,0.0046,0.0054,0.0083,0.0169,0.0097,0.0040,0.0074,0.0041,0.0046,0.0070,0.0112,0.0157,0.0204,0.0185,0.0093,0.0038,0.0101,0.0060,0.0099,0.0180,0.0276,0.0329,0.0321,0.0187,0.0090},
        {0.0045,0.0041,0.0092,0.0044,0.0043,0.0048,0.0057,0.0088,0.0178,0.0102,0.0043,0.0033,0.0092,0.0049,0.0074,0.0118,0.0165,0.0215,0.0195,0.0098,0.0040,0.0045,0.0135,0.0105,0.0190,0.0291,0.0347,0.0339,0.0197,0.0095},
        {0.0046,0.0041,0.0044,0.0096,0.0043,0.0049,0.0058,0.0089,0.0182,0.0104,0.0043,0.0034,0.0044,0.0107,0.0075,0.0120,0.0168,0.0219,0.0199,0.0100,0.0041,0.0046,0.0064,0.0229,0.0194,0.0296,0.0354,0.0345,0.0201,0.0097},
        {0.0045,0.0040,0.0043,0.0043,0.0097,0.0048,0.0057,0.0087,0.0177,0.0101,0.0042,0.0033,0.0042,0.0048,0.0168,0.0117,0.0164,0.0213,0.0193,0.0097,0.0040,0.0045,0.0062,0.0104,0.0435,0.0288,0.0344,0.0336,0.0196,0.0094},
        {0.0051,0.0046,0.0048,0.0049,0.0048,0.0135,0.0065,0.0099,0.0201,0.0115,0.0048,0.0037,0.0048,0.0055,0.0083,0.0329,0.0186,0.0242,0.0220,0.0110,0.0045,0.0051,0.0071,0.0118,0.0214,0.0811,0.0391,0.0381,0.0222,0.0107},
        {0.0060,0.0054,0.0057,0.0058,0.0057,0.0065,0.0199,0.0117,0.0239,0.0136,0.0057,0.0044,0.0057,0.0065,0.0099,0.0158,0.0573,0.0287,0.0261,0.0131,0.0054,0.0061,0.0084,0.0140,0.0255,0.0389,0.1204,0.0453,0.0264,0.0127},
        {0.0092,0.0083,0.0088,0.0089,0.0087,0.0099,0.0117,0.0464,0.0364,0.0208,0.0087,0.0068,0.0087,0.0100,0.0150,0.0241,0.0337,0.1137,0.0398,0.0200,0.0082,0.0093,0.0129,0.0214,0.0388,0.0594,0.0709,0.1793,0.0403,0.0194},
        {0.0187,0.0169,0.0178,0.0182,0.0177,0.0201,0.0239,0.0364,0.1749,0.0424,0.0177,0.0138,0.0178,0.0203,0.0306,0.0491,0.0686,0.0892,0.1914,0.0407,0.0167,0.0188,0.0262,0.0435,0.0790,0.1208,0.1442,0.1407,0.1937,0.0395},
        {0.0107,0.0097,0.0102,0.0104,0.0101,0.0115,0.0136,0.0208,0.0424,0.1355,0.0101,0.0079,0.0102,0.0116,0.0175,0.0281,0.0393,0.0510,0.0464,0.1301,0.0096,0.0108,0.0150,0.0249,0.0452,0.0691,0.0825,0.0805,0.0469,0.1263},
        {0.0142,0.0040,0.0043,0.0043,0.0042,0.0048,0.0057,0.0087,0.0177,0.0101,0.0134,0.0033,0.0042,0.0048,0.0073,0.0117,0.0164,0.0213,0.0194,0.0097,0.0127,0.0045,0.0062,0.0104,0.0189,0.0289,0.0344,0.0336,0.0196,0.0094},
        {0.0035,0.0074,0.0033,0.0034,0.0033,0.0037,0.0044,0.0068,0.0138,0.0079,0.0033,0.0061,0.0033,0.0038,0.0057,0.0092,0.0128,0.0166,0.0151,0.0076,0.0031,0.0083,0.0049,0.0081,0.0147,0.0225,0.0269,0.0262,0.0153,0.0074},
        {0.0045,0.0041,0.0092,0.0044,0.0042,0.0048,0.0057,0.0087,0.0178,0.0102,0.0042,0.0033,0.0092,0.0049,0.0073,0.0118,0.0165,0.0214,0.0195,0.0098,0.0040,0.0045,0.0135,0.0104,0.0190,0.0290,0.0346,0.0338,0.0197,0.0095},
        {0.0051,0.0046,0.0049,0.0107,0.0048,0.0055,0.0065,0.0100,0.0203,0.0116,0.0048,0.0038,0.0049,0.0119,0.0084,0.0134,0.0188,0.0244,0.0222,0.0111,0.0046,0.0052,0.0072,0.0256,0.0216,0.0331,0.0395,0.0385,0.0225,0.0108},
        {0.0077,0.0070,0.0074,0.0075,0.0168,0.0083,0.0099,0.0150,0.0306,0.0175,0.0073,0.0057,0.0073,0.0084,0.0291,0.0203,0.0284,0.0368,0.0335,0.0168,0.0069,0.0078,0.0108,0.0180,0.0753,0.0499,0.0596,0.0581,0.0339,0.0163},
        {0.0124,0.0112,0.0118,0.0120,0.0117,0.0329,0.0158,0.0241,0.0491,0.0281,0.0117,0.0092,0.0118,0.0134,0.0203,0.0806,0.0455,0.0592,0.0538,0.0270,0.0111,0.0125,0.0173,0.0288,0.0524,0.1983,0.0956,0.0933,0.0544,0.0262},
        {0.0173,0.0157,0.0165,0.0168,0.0164,0.0186,0.0573,0.0337,0.0686,0.0393,0.0164,0.0128,0.0165,0.0188,0.0284,0.0455,0.1650,0.0827,0.0751,0.0377,0.0155,0.0175,0.0242,0.0403,0.0733,0.1120,0.3466,0.1304,0.0761,0.0366},
        {0.0225,0.0204,0.0215,0.0219,0.0213,0.0242,0.0287,0.1137,0.0892,0.0510,0.0213,0.0166,0.0214,0.0244,0.0368,0.0592,0.0827,0.2786,0.0976,0.0490,0.0201,0.0227,0.0315,0.0524,0.0952,0.1456,0.1737,0.4394,0.0988,0.0476},
        {0.0204,0.0185,0.0195,0.0199,0.0193,0.0220,0.0261,0.0398,0.1914,0.0464,0.0194,0.0151,0.0195,0.0222,0.0335,0.0538,0.0751,0.0976,0.2096,0.0445,0.0183,0.0206,0.0286,0.0476,0.0865,0.1323,0.1579,0.1540,0.2121,0.0432},
        {0.0103,0.0093,0.0098,0.0100,0.0097,0.0110,0.0131,0.0200,0.0407,0.1301,0.0097,0.0076,0.0098,0.0111,0.0168,0.0270,0.0377,0.0490,0.0445,0.1250,0.0092,0.0104,0.0144,0.0239,0.0434,0.0664,0.0792,0.0773,0.0451,0.1213},
        {0.0134,0.0038,0.0040,0.0041,0.0040,0.0045,0.0054,0.0082,0.0167,0.0096,0.0127,0.0031,0.0040,0.0046,0.0069,0.0111,0.0155,0.0201,0.0183,0.0092,0.0120,0.0043,0.0059,0.0098,0.0178,0.0273,0.0326,0.0318,0.0185,0.0089},
        {0.0047,0.0101,0.0045,0.0046,0.0045,0.0051,0.0061,0.0093,0.0188,0.0108,0.0045,0.0083,0.0045,0.0052,0.0078,0.0125,0.0175,0.0227,0.0206,0.0104,0.0043,0.0113,0.0067,0.0111,0.0201,0.0307,0.0367,0.0358,0.0209,0.0100},
        {0.0066,0.0060,0.0135,0.0064,0.0062,0.0071,0.0084,0.0129,0.0262,0.0150,0.0062,0.0049,0.0135,0.0072,0.0108,0.0173,0.0242,0.0315,0.0286,0.0144,0.0059,0.0067,0.0198,0.0154,0.0279,0.0427,0.0509,0.0497,0.0290,0.0139},
        {0.0110,0.0099,0.0105,0.0229,0.0104,0.0118,0.0140,0.0214,0.0435,0.0249,0.0104,0.0081,0.0104,0.0256,0.0180,0.0288,0.0403,0.0524,0.0476,0.0239,0.0098,0.0111,0.0154,0.0549,0.0464,0.0710,0.0847,0.0826,0.0482,0.0232},
        {0.0199,0.0180,0.0190,0.0194,0.0435,0.0214,0.0255,0.0388,0.0790,0.0452,0.0189,0.0147,0.0190,0.0216,0.0753,0.0524,0.0733,0.0952,0.0865,0.0434,0.0178,0.0201,0.0279,0.0464,0.1945,0.1290,0.1539,0.1502,0.0876,0.0421},
        {0.0305,0.0276,0.0291,0.0296,0.0288,0.0811,0.0389,0.0594,0.1208,0.0691,0.0289,0.0225,0.0290,0.0331,0.0499,0.1983,0.1120,0.1456,0.1323,0.0664,0.0273,0.0307,0.0427,0.0710,0.1290,0.4880,0.2353,0.2296,0.1339,0.0644},
        {0.0364,0.0329,0.0347,0.0354,0.0344,0.0391,0.1204,0.0709,0.1442,0.0825,0.0344,0.0269,0.0346,0.0395,0.0596,0.0956,0.3466,0.1737,0.1579,0.0792,0.0326,0.0367,0.0509,0.0847,0.1539,0.2353,0.7283,0.2740,0.1598,0.0769},
        {0.0355,0.0321,0.0339,0.0345,0.0336,0.0381,0.0453,0.1793,0.1407,0.0805,0.0336,0.0262,0.0338,0.0385,0.0581,0.0933,0.1304,0.4394,0.1540,0.0773,0.0318,0.0358,0.0497,0.0826,0.1502,0.2296,0.2740,0.6933,0.1559,0.0750},
        {0.0207,0.0187,0.0197,0.0201,0.0196,0.0222,0.0264,0.0403,0.1937,0.0469,0.0196,0.0153,0.0197,0.0225,0.0339,0.0544,0.0761,0.0988,0.2121,0.0451,0.0185,0.0209,0.0290,0.0482,0.0876,0.1339,0.1598,0.1559,0.2147,0.0437},
        {0.0100,0.0090,0.0095,0.0097,0.0094,0.0107,0.0127,0.0194,0.0395,0.1263,0.0094,0.0074,0.0095,0.0108,0.0163,0.0262,0.0366,0.0476,0.0432,0.1213,0.0089,0.0100,0.0139,0.0232,0.0421,0.0644,0.0769,0.0750,0.0437,0.1177}};

	// Statistical Errors
    double sigmaStatRemu[] = {0.00077,0.00063,0.00062,0.00053,0.00052,0.00057,0.00061,0.00078,0.00085,0.00103,0.00050,0.00049,0.00046,0.00047,0.00070,0.00113,0.00152,0.00136,0.00077,0.00089,0.00056,
                    0.00056,0.00073,0.00106,0.00180,0.00272,0.00287,0.00221,0.00093,0.00111};

    // Put the errors together in a matrix
    TMatrix sigmaRemu_inv(30,30);
    TMatrix sigmaRemu(30,30);
    for(int i = 0; i < 30; i ++){
        for(int j = 0; j < 30; j++){
            if(i == j)
                sigmaRemu_inv(i,j) = sigmaSystRemu[i][j] + pow(sigmaStatRemu[i],2) * pow(10,5);
            else
                sigmaRemu_inv(i,j) = sigmaSystRemu[i][j];
        }
    }
    sigmaRemu = sigmaRemu_inv.Invert();
    for(int i = 0; i < 30; i++){
        for(int j = 0; j < 30; j++){
            sigmaRemu[i][j] *= pow(10,5);
            pack.sigmaRemu[i][j] = sigmaRemu[i][j];
        }
    }

    double energyMin[] = {3.,12.,16.,20.,25.,30.,35.,40.,50.,100.};
    double energyMax[] = {12.,16.,20.,25.,30.,35.,40.,50.,100.,170.};

    double fullMixHighDm2Remu[] = {0.4406,0.4327,0.4055,0.4073,0.3979,0.3968,0.4030,0.3823,0.3636,0.3427,0.4252,0.3991,0.4041,0.3914,0.4028,0.3887,0.3771,0.3777,0.3605,0.3186,0.4292,0.4187,0.4138,0.4075,
                    0.3976,0.4013,0.4003,0.3795,0.3797,0.2960};

    double l1 = .422; double l2 = .837;
    pack.norm = new double[30];

    // First, find the norm
    for(int i = 0; i < 30; i++)
        pack.norm[i] = 0;

    for(int i = 0; i < maxEnergyBins; i++){
        double EMin = energyMin[i];
        double EMax = energyMax[i];
        for(int j = 0; j < maxRadialBins; j++){
            int iN = maxEnergyBins*j + i;
            double res = (l2 - l1) * (EMax - EMin);
            pack.norm[iN] = 2. * fullMixHighDm2Remu[iN] / res;
        }
    }

    // Now, fill the sinsq grid vectors
    for(int i = 0; i < maxEnergyBins; i++){
        double EMin = energyMin[i];
        double EMax = energyMax[i];
        integralFuncs._EnuAvg = (energyMax[i] + energyMin[i])/2.;
        for(int j = 0; j < maxRadialBins; j++){
            int iN = maxEnergyBins*j + i;
            for(int k = 0; k < dm2VecMaxDim; k++){
                integralFuncs._dm2 = dm2Vec[k];

                ROOT::Math::Functor1D wf1(&integralFuncs, &integralFuncsNomad::sinSqFunction);
                ROOT::Math::Functor1D wf2(&integralFuncs, &integralFuncsNomad::sinSqFunctionCPV);
                ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
                ig.SetAbsTolerance(1*pow(10,-4));
                ig.SetRelTolerance(1*pow(10,-4));

                ig.SetFunction(wf1);
                pack.sinSqDeltaGrid[k][iN] = ig.Integral(EMin, EMax);// (EMax - EMin);
                ig.SetFunction(wf2);
                pack.sinSqDeltaGrid2[k][iN] = ig.Integral(EMin, EMax);// (EMin - EMax);
			}
        }
    }

    ndf += maxEnergyBins * maxRadialBins;
	if(debug) std::cout << "Nomad initialized. Bins: " << maxEnergyBins * maxRadialBins << std::endl;
    return pack;
    }
ccfrPackage ccfrInit(){
	ccfrPackage pack;

	const int maxEnergyBins = 18;

	double temp0[] = {0.966,1.014,1.039,0.949,0.988,1.026,1.017,0.963,0.959,1.018,0.990,1.021,1.027,0.971,1.028,1.019,0.938,0.945};
	pack.observed = new double[maxEnergyBins];
	for(int i = 0; i < maxEnergyBins; i++)
		pack.observed[i] = temp0[i];

	pack.m_front = 105.; 	pack.m_back = 444.;

	pack.sigmaRatio.resize(maxEnergyBins, std::vector<double>(maxEnergyBins));
	pack.sinSqDeltaGrid_front.resize(dm2VecMaxDim, std::vector<double>(maxEnergyBins));
	pack.sinSqDeltaGrid_back.resize(dm2VecMaxDim, std::vector<double>(maxEnergyBins));
	pack.noOscGrid.resize(maxEnergyBins, std::vector<double>(2));

	TMatrix sigmaRatio(maxEnergyBins, maxEnergyBins);
	TMatrix sigmaRatio_inv(maxEnergyBins,maxEnergyBins);
	double NuEnergy[] = {40.76,36.42,96.08,53.62,45.20,128.50,61.01,49.78,151.28,70.51,55.05,176.25,82.86,63.00,209.60,61.06,49.46,150.18};
    double sigmaRatio1[] = {0.034,0.032,0.057,0.025,0.026,0.036,0.032,0.033,0.036,0.025,0.029,0.027,0.030,0.037,0.028,0.032,0.031,0.060};
	double sigmaSyst = 0.015;
	for(int i = 0; i < maxEnergyBins; i++){
    	for(int j = 0; j < maxEnergyBins; j++){
        	if(i == j)  sigmaRatio_inv[i][j] = pow(sigmaRatio1[i],2) + pow(sigmaSyst,2);
            else    sigmaRatio_inv[i][j] = pow(sigmaSyst,2);
		}
	}

	sigmaRatio = sigmaRatio_inv.Invert();
	for(int i = 0; i < maxEnergyBins; i++){
		for(int j = 0; j < maxEnergyBins; j++){
			pack.sigmaRatio[i][j] = sigmaRatio[i][j];
		}
	}

	double lBack = 1.116; double lFront = 0.715;
	double deltaLBack = .357; double deltaLFront = .3565;
	double l1Back = lBack - deltaLBack/2;   double l2Back = lBack + deltaLBack/2;
	double l1Front = lFront - deltaLFront/2;double l2Front = lFront + deltaLFront/2;

	// Front Detector
	for(int iC = 0; iC < maxEnergyBins; iC++){
	    pack.noOscGrid[iC][0] = (l2Front - l1Front) / (l1Front * l2Front);
	}
	for(int k = 0; k < dm2VecMaxDim; k ++){
		for(int iC = 0; iC < maxEnergyBins; iC++){
		    double Enu = NuEnergy[iC];
		    double delta = 1.27 * dm2Vec[k] / Enu;
		    double num1 = 1. - cos(2. * delta * l1Front) - 2. * delta * l1Front * sineInt(2. * delta * l1Front);
		    double num2 = 1. - cos(2. * delta * l2Front) - 2. * delta * l2Front * sineInt(2. * delta * l2Front);
		    double den1 = 2 * l1Front;  double den2 = 2* l2Front;
		    double term1 = num1 / den1; double term2 = num2 / den2;
		    pack.sinSqDeltaGrid_front[k][iC] = term1 - term2;
		}
	}
	// Back Detector
	for(int iC = 0; iC < maxEnergyBins; iC++){
	    pack.noOscGrid[iC][1] = (l2Back - l1Back) / (l1Back * l2Back);
	}
	for(int k = 0; k < dm2VecMaxDim; k ++){
	    for(int iC = 0; iC < maxEnergyBins; iC++){
	        double Enu = NuEnergy[iC];
	        double delta = 1.27 * dm2Vec[k] / Enu;
	        double num1 = 1. - cos(2. * delta * l1Back) - 2. * delta * l1Back * sineInt(2. * delta * l1Back);
			double num2 = 1. - cos(2. * delta * l2Back) - 2. * delta * l2Back * sineInt(2. * delta * l2Back);
	        double den1 = 2 * l1Back;  double den2 = 2* l2Back;
	        double term1 = num1 / den1; double term2 = num2 / den2;
	        pack.sinSqDeltaGrid_back[k][iC] = term1 - term2;
		}
	}

	ndf += maxEnergyBins;
	if(debug) std::cout << "CCFR initialized. Bins: " << maxEnergyBins << std::endl;
    return pack;
}

double getMuEnergy(double range){
	// For CDHS, we have muon ranges in iron and we need to convert them to energy
	double betaGammaVec[] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
	double rangeOverMassVec[] = {100.,600.,1100.,2000.,2500.,3100.,4000.,4500.,5000.,6000.,11000,18000.,22000.,28000.,32000.,39000.,43000,500000.,530000.};
    double ironDensity = 7.87;
	double muMass = .105;

	// c    there is a difference between the projected muon range measured by CDHS
	// c    and the full range first, the two detectors were at at 22deg angle wrt
	// c    nu beamline; second, a muon typically comes off at a 20deg angle from a
	// c    neutrino interaction
	double scattAngle = 20. * TMath::Pi()/180.;
	double detAngle = 22. * TMath::Pi()/180.;

	double project = cos(scattAngle) * cos(detAngle);
	double rangeOverMass = range * ironDensity / muMass;
	double betaGamma;

	for(int i = 1; i < 19; i ++){
		if(rangeOverMass >= rangeOverMassVec[i-1] * project && rangeOverMass < rangeOverMassVec[i] * project){
			betaGamma = betaGammaVec[i-1] + (betaGammaVec[i] - betaGammaVec[i-1]) * (rangeOverMass - rangeOverMassVec[i-1]) / (rangeOverMassVec[i] - rangeOverMassVec[i-1]);
			break;
		}
	}
    if(rangeOverMass < rangeOverMassVec[0] * project)
		betaGamma = betaGammaVec[0];
	else if (rangeOverMass >= rangeOverMassVec[18] * project)
		betaGamma = betaGammaVec[18];
	double gamma = sqrt(1. + pow(betaGamma,2));

	return gamma * muMass;
}
cdhsPackage cdhsInit(){
	cdhsPackage pack;

  	const int nBins = 15;

  	pack.sigmaRatio.resize(nBins, std::vector<double>(nBins));
	pack.dm2Vec.resize(601);
	pack.sinSqDeltaGrid_front.resize(dm2VecMaxDim, std::vector<double>(nBins));
	pack.sinSqDeltaGrid_back.resize(dm2VecMaxDim, std::vector<double>(nBins));
	pack.noOscGrid.resize(nBins, std::vector<double>(2));

	double temp0[] = {0.985,1.006,0.968,1.148,1.000,1.137,1.155,0.887,1.123,0.973,1.039,1.187,1.196,1.006,0.963};
	double temp1[] = {750.,750.,750.,750.,750.,750.,750.,750.,250.,250.,250.,250.,250.,250.,250.};
	double temp2[] = {150.,150.,150.,150.,150.,150.,150.,150.,100.,100.,100.,100.,100.,100.,100.};
   	pack.observed = new double[nBins];
   	pack.m_front = new double[nBins];
   	pack.m_back = new double[nBins];
   	for(int i = 0; i < nBins; i++){
   		pack.observed[i] = temp0[i];
		pack.m_front[i] = temp2[i];
		pack.m_back[i] = temp1[i];
   	}

    double sigmaRatio1[] = {0.066,0.055,0.070,0.075,0.087,0.104,0.135,0.200,0.078,0.075,0.085,0.103,0.120,0.132,0.110};

    TMatrix sigmaRatio(nBins,nBins);
    TMatrix sigmaRatio_inv(nBins,nBins);
    double sigmaSyst = 0.025;
    for(int i = 0; i < nBins; i++){
        for(int j = 0; j < nBins; j++){
            if(i == j)  sigmaRatio_inv[i][j] = pow(sigmaRatio1[i],2) + pow(sigmaSyst,2);
            else    sigmaRatio_inv[i][j] = pow(sigmaSyst,2);
        }
    }
    sigmaRatio = sigmaRatio_inv.Invert();
    for(int i = 0; i < nBins; i++){
    	for(int j = 0; j < nBins; j++){
    		pack.sigmaRatio[i][j] = sigmaRatio[i][j];
    	}
    }

	// The cdhs integrals take forever, so i just ran them from the original fortran code and stored them in a txt file, which i read here!
	ifstream file;
	file.open(dataLoc+"cdhs_sinsq.txt");
    for(int i = 0; i < nBins; i++){
		for(int k = 0; k < 601; k++){
			file >> pack.dm2Vec[k];
			file >> pack.sinSqDeltaGrid_front[k][i];
			file >> pack.sinSqDeltaGrid_back[k][i];
			//std::cout << pack.dm2Vec[k] << " " << pack.sinSqDeltaGrid_front[k][i] << " " << pack.sinSqDeltaGrid_back[k][i] << std::endl;
		}
    }
    file.close();
	file.open(dataLoc+"cdhs_noosc.txt");
    for(int i = 0; i < nBins; i++){
		file >> pack.noOscGrid[i][0];
		file >> pack.noOscGrid[i][1];
    }
    file.close();

	/*
    double minRange[15] = {37.5,50.,75.,100.,137.5,187.5,250.,337.5,40.,50.,70.,100.,125.,170.,225.};
    double maxRange[15] = {50.,75.,100.,137.5,187.5,250.,337.5,450.,50.,70.,100.,125.,170.,225.,450.};
    double lBack = .885; double lFront = .130;
    double deltaLBack = .144; double deltaLFront = .074;

    double l1Back = lBack - deltaLBack/2;   double l2Back = lBack + deltaLBack/2;
    double l1Front = lFront - deltaLFront/2;double l2Front = lFront + deltaLFront/2;
    integralFuncsCDHS2 integralFuncs;

    for(int i = 0; i < nBins; i++){
        integralFuncs._noOsc = true;    integralFuncs._front = true;
        ROOT::Math::Functor1D wf(&integralFuncs, &integralFuncsCDHS2::intFunc2);
        ROOT::Math::Integrator ig;
        ig.SetFunction(wf);
		double emuMin = getMuEnergy(minRange[i]);
		double emuMax = getMuEnergy(maxRange[i]);
        pack.noOscGrid[i][0] = ig.Integral(emuMin,emuMax);

        for(int k = 0; k < dm2VecMaxDim; k++){
			integralFuncs._noOsc = false;    integralFuncs._dm2 = dm2Vec[k];
            pack.sinSqDeltaGrid_front[k][i] = ig.Integral(emuMin, emuMax);// / (emuMax - emuMin);
			//std::cout << "sinsq front: " << pack.sinSqDeltaGrid_front[k][i] << std::endl;
		}
    }
    for(int i = 0; i < nBins; i++){
        integralFuncs._noOsc = true;    integralFuncs._front = false;
        ROOT::Math::Functor1D wf(&integralFuncs, &integralFuncsCDHS2::intFunc2);
        ROOT::Math::Integrator ig;
        ig.SetFunction(wf);
		double emuMin = getMuEnergy(minRange[i]);
		double emuMax = getMuEnergy(maxRange[i]);
        pack.noOscGrid[i][1] = ig.Integral(emuMin,emuMax);

        for(int k = 0; k < dm2VecMaxDim; k++){
            integralFuncs._noOsc = false;    integralFuncs._dm2 = dm2Vec[k];
            pack.sinSqDeltaGrid_back[k][i] = ig.Integral(emuMin, emuMax);// / (emuMax - emuMin);
			//std::cout << pack.sinSqDeltaGrid_back[k][i] << std::endl;
		}
    }
	*/

    ndf += nBins;
	if(debug) std::cout << "CDHS initialized. Bins: " << nBins << std::endl;
    return pack;
}
bugeyPackage bugeyInit(){
  	bugeyPackage pack;

  	int maxEnergyBins = 25;
  	int nBaselines = 3;

  	double deltaE[] = {.2, .2, .5};
  	double EMin = 1.; double EMax = 6.;

  	pack.sigmaBigA = 4.796e-2;
  	pack.sigmaB =2.e-2;
  	pack.sigmaSmallA =1.414e-2;

  	pack.observed.resize(nBaselines, std::vector<double>(maxEnergyBins));
  	pack.sigmaRatio.resize(nBaselines, std::vector<double>(maxEnergyBins));
  	pack.energy.resize(nBaselines, std::vector<double>(maxEnergyBins));
  	pack.sinSqDeltaGrid.resize(dm2VecMaxDim, std::vector<std::vector<double>>(maxEnergyBins, std::vector<double>(nBaselines)));

  	double ratio_obs1[] = {1.0182,1.0017,0.9815,0.9875,0.9981,0.9849,0.9749,0.9743,0.9840,0.9910,1.0042,0.9966,0.9629,1.0144,
      	0.9971,0.9759,0.9889,0.9493,1.0056,0.8928,1.0160,0.9291,0.8585,1.0012,0.8981};
  	double ratio_obs2[] = {0.9589,1.0344,0.9899,1.0071,0.9990,1.0267,0.9830,1.0298,0.9889,0.9441,0.9827,1.0123,1.0069,0.9059,
      	0.9877,1.0281,1.0104,0.8905,1.0377,0.9940,1.0593,0.9457,0.9146,0.9941,0.9161};
  	double ratio_obs3[] = {0.1685,0.8189,1.2169,1.3342,0.7995,1.1623,1.2918,1.3251,1.2711,0.5812};

  	double sigmaRatio1[] = {0.01893,0.01610,0.01516,0.01472,0.01461,0.01463,0.01522,0.01595,0.01660,0.01906,0.02061,0.02064,
      	0.02071,0.02366,0.02498,0.02681,0.02914,0.03117,0.03657,0.03699,0.04789,0.05124,0.05401,0.07408,0.07716};
  	double sigmaRatio2[] = {0.06553,0.04876,0.03628,0.03182,0.03129,0.03140,0.03128,0.03082,0.02922,0.02847,0.03159,0.03412,
      	0.03705,0.03669,0.04426,0.04753,0.05371,0.05415,0.06137,0.07084,0.07786,0.08788,0.10928,0.13213,0.17449};
  	double sigmaRatio3[] = {0.54296,0.29210,0.20995,0.18166,0.23880,0.26919,0.34341,0.28928,0.79799,0.90417};

  	double normReactorAno[] = {1.06237,1.06197,1.0627};
  	double nBins[] = {25, 25, 10};

  	double baselines[] = {15., 40., 95.};

  	for(int j = 0; j < nBaselines; j++){
    	for(int i = 0; i < nBins[j]; i++){
      		pack.energy[j][i] = EMin + (double(i) + 0.5)/(nBins[j])*(EMax - EMin);
      		if(j==0){
        		pack.observed[j][i] = ratio_obs1[i];
        		pack.sigmaRatio[j][i] = sigmaRatio1[i];
      		}
      		if(j==1){
        		pack.observed[j][i] = ratio_obs2[i];
        		pack.sigmaRatio[j][i] = sigmaRatio2[i];
      		}
      		if(j==2){
        		pack.observed[j][i] = ratio_obs3[i];
        		pack.sigmaRatio[j][i] = sigmaRatio3[i];
      		}
    	}
  	}

  	// Now, for the more complicated stuff
  	double EnuMin = 1.8; double EnuMax = 10.;

  	integralFuncsBugey integralFuncs;
  	ROOT::Math::Functor1D wf(&integralFuncs, &integralFuncsBugey::sinSqFunction);
  	ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  	ig.SetAbsTolerance(1*pow(10,-4));
  	ig.SetRelTolerance(1*pow(10,-4));
  	ig.SetFunction(wf);

  	for(int k = 0; k < dm2VecMaxDim; k++){
    	for(int j = 0; j < nBaselines; j++){
      		for(int i = 0; i < nBins[j]; i++){

        		double Enu = pack.energy[j][i] + 1.8;
        		int n = 1.27 * dm2Vec[k] * baselines[j] / (Enu * TMath::Pi());
        		double eNuNext = 1.27 * dm2Vec[k] * baselines[j] / ((n + 0.5) * TMath::Pi());

        		if(abs(Enu - eNuNext) >= 0.123 * sqrt(Enu - 1.8)){
          			integralFuncs._energy = Enu;
          			integralFuncs._dm2 = dm2Vec[k];
          			integralFuncs._jB = j;

          			pack.sinSqDeltaGrid[k][i][j] = ig.Integral(EnuMin,EnuMax);// / (EnuMax-EnuMin);
        		}
        		else
          			pack.sinSqDeltaGrid[k][i][j] = .5;
      		}
    	}
  	}
	ndf += 60;
	if(debug) std::cout << "Bugey initialized. Bins: " << 60 << std::endl;
  	return pack;
}
choozPackage choozInit(){

	choozPackage pack;

	const int maxEnergyBins = 7;
	const int nBaselines = 2;

	double EMin = 0.8; double EMax = 6.4;
	double baselines[2] = {910.,1025.};

	// Positron yield in a given (baseline,energy) bin
	double obs1[7] = {0.151,0.490,0.656,0.515,0.412,0.248,0.102};
	double obs2[7] = {0.176,0.510,0.610,0.528,0.408,0.231,0.085};

	// error matrix from statistical error on positron yield
	double sigmax1[7] = {0.031,0.039,0.041,0.036,0.033,0.030,0.023};
	double sigmax2[7] = {0.035,0.047,0.049,0.044,0.040,0.034,0.026};
	double sigmax12[7] = {-2.2e-4,-1.5e-4,-3.5e-4,-3.3e-4,-2.0e-4,-0.7e-4,-1.3e-4};

	double reactorAno = 1.051;

	pack.deltae = 0.8;
	pack.sigmaAlpha = 2.7e-2;
	pack.sigmaG = 1.1e-2;

	if(ReactorAnomaly)
		pack.sigmaAlpha = 3.3e-2;

	pack.energy.resize(7);
	pack.observed.resize(14);

	for(int i = 0; i < 7; i++){
		pack.energy[i] = EMin + (double(i)+0.5)/maxEnergyBins;
		pack.observed[i] = obs1[i];
		pack.observed[i+maxEnergyBins] = obs2[i];
	}

	double _noOsc[14] = {0.172,0.532,0.632,0.530,0.379,0.208,0.101,0.172,0.532,0.632,0.530,0.379,0.208,0.101};
	pack.noOsc.resize(14);
	for(int i = 0; i < 14; i++)
		pack.noOsc[i] = _noOsc[i];

	if(ReactorAnomaly){
		for(int j = 0; j < 14; j++){
			pack.noOsc[j] *= reactorAno;
		}
	}

	std::vector <double> sigmaSyst;
		sigmaSyst.push_back(4.0e-2*pack.noOsc[0]);
		sigmaSyst.push_back(2.0e-2*pack.noOsc[1]);
		sigmaSyst.push_back(2.0e-2*pack.noOsc[2]);
		sigmaSyst.push_back(2.5e-2*pack.noOsc[3]);
		sigmaSyst.push_back(4.0e-2*pack.noOsc[4]);
		sigmaSyst.push_back(5.5e-2*pack.noOsc[5]);
		sigmaSyst.push_back(12.0e-2*pack.noOsc[6]);
		sigmaSyst.push_back(4.0e-2*pack.noOsc[7]);
		sigmaSyst.push_back(2.0e-2*pack.noOsc[8]);
		sigmaSyst.push_back(2.0e-2*pack.noOsc[9]);
		sigmaSyst.push_back(2.5e-2*pack.noOsc[10]);
		sigmaSyst.push_back(4.0e-2*pack.noOsc[11]);
		sigmaSyst.push_back(5.5e-2*pack.noOsc[12]);
		sigmaSyst.push_back(12.0e-2*pack.noOsc[13]);

	// Define error matrix for positron yields
	TMatrix sigmaXM(14,14);
	pack.sigmaX.resize(14, std::vector<double>(14));
	for(int i = 0; i < 14; i++){
		for(int j = 0; j < 14; j++){
			if(i == j){
				if(i < 7)
					sigmaXM[i][j] = pow(sigmax1[i],2) + pow(sigmaSyst[i],2);
				else
					sigmaXM[i][j] = pow(sigmax2[i-maxEnergyBins],2) + pow(sigmaSyst[i],2);
			}
			else if(i == j-7 || i == j+7){
				if(i < 7)
					sigmaXM[i][j] = sigmax12[i];
				else
					sigmaXM[i][j] = sigmax12[j];
			}
			else
				sigmaXM[i][j] = 0;
		}
	}
	sigmaXM.Invert();
	for(int i = 0; i < 14; i ++)
		for(int j = 0; j < 14; j++)
			pack.sigmaX[i][j] = sigmaXM[i][j];

	return pack;
}
xsecPackage xsecInit(){

	xsecPackage pack;

	pack.karmen_Enu.push_back(28.70);
	pack.karmen_Enu.push_back(32.71);
	pack.karmen_Enu.push_back(36.55);
	pack.karmen_Enu.push_back(40.90);
	pack.karmen_Enu.push_back(45.42);
	pack.karmen_Enu.push_back(49.51);

	pack.karmen.push_back(4.68);
	pack.karmen.push_back(6.31);
	pack.karmen.push_back(9.18);
	pack.karmen.push_back(13.17);
	pack.karmen.push_back(23.57);
	pack.karmen.push_back(36.77);

	pack.karmen_error.push_back(1.26);
	pack.karmen_error.push_back(1.39);
	pack.karmen_error.push_back(1.52);
	pack.karmen_error.push_back(1.83);
	pack.karmen_error.push_back(2.47);
	pack.karmen_error.push_back(5.76);

	pack.lsnd_Enu.push_back(37.80);
	pack.lsnd_Enu.push_back(41.48);
	pack.lsnd_Enu.push_back(44.93);
	pack.lsnd_Enu.push_back(47.89);
	pack.lsnd_Enu.push_back(49.91);

	pack.lsnd.push_back(11.42);
	pack.lsnd.push_back(17.59);
	pack.lsnd.push_back(21.39);
	pack.lsnd.push_back(25.20);
	pack.lsnd.push_back(35.17);

	pack.lsnd_error.push_back(0.92);
	pack.lsnd_error.push_back(1.31);
	pack.lsnd_error.push_back(1.58);
	pack.lsnd_error.push_back(2.23);
	pack.lsnd_error.push_back(4.99);

	pack.ESigma = 0.15;

	pack.lsnd_sys = sqrt(pow(.1,2) - pow(.07,2));
	pack.karmen_sys = sqrt(pow(.088,2) - pow(.07,2));
	pack.correl_sys = sqrt(pow(.12,2) + pow(.07,2));


	ndf += 11;
	if(debug) std::cout << "XSEC initialized. Bins: " << 11 << std::endl;
	return pack;
}
