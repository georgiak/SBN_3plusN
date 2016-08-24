/* ------------------------------------------//
Created by Davio Cianci and Georgia Karagiorgi
Jan 25th, 2016

Notes:
 - found error in fortran code in ourpred_ws for minos

------------------------------------------// */

#include "globalFit.h"
#include <time.h>

// Initializations
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "detectorInit.cpp+"
#else
#include "detectorInit.cpp"
#endif

// ChiSq Calculations
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "detectorChisq.cpp+"
#else
#include "detectorChisq.cpp"
#endif

// Declare Mass and Mixing Variables
float ran[13];
int iMCGen;             // Index for MC model
double dm2Min[3], dm2Max[3];
double dm2Minimum, dm2Maximum;
float UMin = 0;

bool reject1, reject2, reject3, reject4, reject;
double chi2Min, chi2LogMin;

int noOfParameters;

TRandom RanGen;

// things for chisq calculation
float m2LogL, m2LogL_det;
float _m2LogL, _m2LogL_det;  // Temporary chisq holders for individual detector calcs

double dm2_, aMuE_, aMuE_CPV_;
double chi2Log, chi2LogOld;

// Ntuple Variables
float chi2, dof, step, temp, gof, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56;

boonePackage mbNuPack, mbNubarPack; atmPackage atmPack; numiPackage numiPack; sinSqPackage lsndPack, karmenPack; galPackage galPack; cdhsPackage cdhsPack;
minosPackage minosPack; minosncPackage minosncPack; booneDisPackage mbNuDisPack, mbNubarDisPack; nomadPackage nomadPack; ccfrPackage ccfrPack;
bugeyPackage bugeyPack; choozPackage choozPack; xsecPackage xsecPack;

bool debug = true;

int globInit(){

    using namespace std;

	jobOptLoc = "/Users/dcianci/Physics/SBN_3plusN/GlobalFits/inputs/"; // /pnfs/lar1nd/scratch/users/dcianci/fits/";
	dataLoc = "/Users/dcianci/Physics/SBN_3plusN/GlobalFits/data/"; ///pnfs/lar1nd/scratch/users/dcianci/fits/data";

    // read jobOption file and fill variables
    jobOpt();
    dm2Max[0] = 100.;   dm2Max[1] = 100.;    dm2Max[2] = 100.;
    dm2Min[0] = .01;     dm2Min[1] = .01;    dm2Min[2] = .01;
	UMax = .5;

    // INITIALIZATIONS
	std::cout << "Start initializations!" << std::endl;

    dm2VecInit(.01, 100.);
	if(MBProcess) mbNuPack = mbNuInit();
	if(debug && MBProcess) std::cout << "MB initialized." << std::endl;
	if(MBProcessNubar) mbNubarPack = mbNubarInit();
	if(debug && MBProcessNubar) std::cout << "MBProce initialized." << std::endl;
	if(ATMOSPHERICProcess) atmPack = atmInit();
	if(debug && ATMOSPHERICProcess) std::cout << "ATMOSPHERIC initialized." << std::endl;
	if(NUMIProcess) numiPack = numiInit();
	if(debug && NUMIProcess) std::cout << "NUMI initialized." << std::endl;
	if(LSNDProcess) lsndPack = lsndInit();
	if(debug && LSNDProcess) std::cout << "LSND initialized." << std::endl;
	if(KARMENProcess) karmenPack = karmenInit();
	if(debug && KARMENProcess) std::cout << "KARMEN initialized." << std::endl;
	if(GALLIUMProcess) galPack = galInit();
	if(debug && GALLIUMProcess) std::cout << "GALLIUM initialized." << std::endl;
	if(MINOSProcess) minosPack = minosInit();
	if(debug && MINOSProcess) std::cout << "MINOS initialized." << std::endl;
	if(MINOSNCProcess) minosncPack = minosncInit();
	if(debug && MINOSNCProcess) std::cout << "MINOSNC initialized." << std::endl;
	if(MBDISProcess) mbNuDisPack = mbNuDisInit();
	if(debug && MBDISProcess) std::cout << "MBDIS initialized." << std::endl;
	if(MBDISProcessNubar) mbNubarDisPack = mbNubarDisInit();
	if(debug && MBDISProcessNubar) std::cout << "MBDISProce initialized." << std::endl;
	if(NOMADProcess) nomadPack = nomadInit();
	if(debug && NOMADProcess) std::cout << "NOMAD initialized." << std::endl;
	if(CCFRProcess) ccfrPack = ccfrInit();
	if(debug && CCFRProcess) std::cout << "CCFR initialized." << std::endl;
	if(CDHSProcess) cdhsPack = cdhsInit();
	if(debug && CDHSProcess) std::cout << "CDHS initialized." << std::endl;
	if(BugeyProcess) bugeyPack = bugeyInit();
	if(debug && BugeyProcess) std::cout << "Bugey initialized." << std::endl;
	if(CHOOZProcess) choozPack = choozInit();
	if(debug && CHOOZProcess) std::cout << "CHOOZ initialized." << std::endl;
	if(XSECProcess) xsecPack = xsecInit();
	if(debug && XSECProcess) std::cout << "XSEC initialized." << std::endl;

    getNDF();
	if(XSECProcess || BugeyProcess || CHOOZProcess)	myMinInit();

	std::cout << "DOF: " << ndf << std::endl;

    std::cout << "Alright! Detector stuff successfully initialized!" << std::endl;

	return 1;
}

int globChisq(int ind){

	std::string outfile = "globPhit.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

	TNtuple *chi2Nt = new TNtuple("chi2Nt","chi2Nt","chi2:step:temp:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

    // Initialize the parameters we'll be using
    neutrinoModel nuModel;
    neutrinoModel nuModelOld;

    // Initialize the chisq result
    chisqStruct chisqDetector, chisqTotal;

	// We're doing a grid scan, motherfuckers
	int nGrid = 100;

	for(int iphi = 0; iphi < nGrid; iphi++){

		nuModel.zero();

		nuModel.Ue[0] = .15; 	nuModel.Um[0] = .17;	nuModel.mNu[0] = sqrt(.92);
		nuModel.Ue[1] = .069;	nuModel.Um[1] = .16; 	nuModel.mNu[1] = sqrt(17);
		nuModel.phi[0] = iphi/float(nGrid) * 2 * TMath::Pi();

        chisqTotal.zero();  chisqDetector.zero();

		clock_t t, st;
  		t = clock();

        // Now, let's actually start calculating the chisq
		if(MINOSProcess == 1){
            chisqDetector = getChi2Minos(nuModel, minosPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Minos: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(CCFRProcess == 1){
			chisqDetector = getChi2CCFR(nuModel, ccfrPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "CCFR: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(CDHSProcess == 1){
			chisqDetector = getChi2CDHS(nuModel, cdhsPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "CDHS: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(MBProcess == 1){
            chisqDetector = getChi2Boone(nuModel, mbNuPack, false);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

			chisqTotal.chi2 += chisqDetector.chi2;
            chisqTotal.chi2_det += chisqDetector.chi2_det;
			if(debug) std::cout << "MB: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(MBProcessNubar == 1){
			chisqDetector = getChi2Boone(nuModel, mbNubarPack, true);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

			chisqTotal.chi2 += chisqDetector.chi2;
			chisqTotal.chi2_det += chisqDetector.chi2_det;
			if(debug) std::cout << "MBNubar: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(LSNDProcess == 1){
			chisqDetector = getLogLikelihood(nuModel, 5, lsndPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "LSND: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(NOMADProcess == 1){
            chisqDetector = getChi2Nomad(nuModel, nomadPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Nomad: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(KARMENProcess == 1){
            chisqDetector = getLogLikelihood(nuModel, 9, karmenPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Karmen: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(NUMIProcess == 1){
            chisqDetector = getChi2Numi(nuModel, numiPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Numi: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(GALLIUMProcess == 1){
            chisqDetector = getChi2Gallium(nuModel, galPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Gal: " << chisqDetector.chi2 << std::endl;
        }



		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(MBDISProcess == 1){
			std::cout << " Here obvi" << std::endl;
			chisqDetector = getChi2MBDis(nuModel, mbNuDisPack);
			std::cout << "chisqdet: " << chisqDetector.chi2 << std::endl;
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
			std::cout << "bannaan" << std::endl;

			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "MBDis: " << chisqDetector.chi2 << std::endl;
		}
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(MBDISProcessNubar == 1){
            chisqDetector = getChi2MBDis(nuModel, mbNubarDisPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "MBDisNubar: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(BugeyProcess == 1){
            chisqDetector = getChi2Bugey(nuModel, bugeyPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Bugey: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(XSECProcess == 1){
            chisqDetector = getChi2Xsec(nuModel, xsecPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Xsec: " << chisqDetector.chi2 << std::endl;
        }
		std::cout << clock() - t << " ticks" << std::endl; t = clock();
		if(chisqTotal.chi2 > chi2Cut) continue;

		if(chisqTotal.chi2 < 0.) {
			std::cout << "WHOAWHOAWHOAAAH HOLD THE GODDAMN PHONE" << std::endl;
			return 0;
		}
        chi2Log = chisqTotal.chi2;

        // Fill Ntuple
        chi2 = chisqTotal.chi2;
        dof = ndf;  m4 = nuModel.mNu[0];    m5 = nuModel.mNu[1];    m6 = nuModel.mNu[2];    ue4 = nuModel.Ue[0];    ue5 = nuModel.Ue[1];    ue6 = nuModel.Ue[2];
        um4 = nuModel.Um[0];    um5 = nuModel.Um[1];    um6 = nuModel.Um[2]; phi45 = nuModel.phi[0];    phi46 = nuModel.phi[1]; phi56 = nuModel.phi[2];
        chi2Nt->Fill(chi2, step, temp, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);
    }

    // Save Ntuple to file
    chi2Nt->Write();
    f->Close();

    return 0;
}


// Since the osc probability is something that'll be played with quite a bit (ie, when we add matter effects) let's put it over here!
oscContribution getOscContributionsNueApp(neutrinoModel model, bool nubar, bool cpv){
	oscContribution oscCon;

	oscCon.dm2[0] = pow(model.mNu[0],2);
	oscCon.aMuE[0] = 4.*model.Ue[0]*model.Um[0]*(model.Ue[0]*model.Um[0] + model.Ue[1]*model.Um[1]*cos(model.phi[0]) + model.Ue[2]*model.Um[2]*cos(model.phi[1]));
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[0] = 2.*model.Ue[0]*model.Um[0]*(model.Ue[1]*model.Um[1]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[1]));
		else    oscCon.aMuE_CPV[0] = -2.*model.Ue[0]*model.Um[0]*(model.Ue[1]*model.Um[1]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[1]));
	}

	oscCon.dm2[1] = pow(model.mNu[1],2);
	oscCon.aMuE[1] = 4.*model.Ue[1]*model.Um[1]*(model.Ue[0]*model.Um[0]*cos(model.phi[0]) + model.Ue[1]*model.Um[1] + model.Ue[2]*model.Um[2]*cos(model.phi[2]));
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[1] = 2.*model.Ue[1]*model.Um[1]*(-model.Ue[0]*model.Um[0]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[2]));
		else    oscCon.aMuE_CPV[1] = -2.*model.Ue[1]*model.Um[1]*(-model.Ue[0]*model.Um[0]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[2]));
	}

	oscCon.dm2[2] = pow(model.mNu[2],2);
	oscCon.aMuE[2] = 4.*model.Ue[2]*model.Um[2]*(model.Ue[0]*model.Um[0]*cos(model.phi[1]) + model.Ue[1]*model.Um[1]*cos(model.phi[2]) + model.Ue[2]*model.Um[2]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[2] = 2.*model.Ue[2]*model.Um[2]*(-model.Ue[0]*model.Um[0]*sin(model.phi[1]) - model.Ue[1]*model.Um[1]*sin(model.phi[2]));
		else    oscCon.aMuE_CPV[2] = -2.*model.Ue[2]*model.Um[2]*(-model.Ue[0]*model.Um[0]*sin(model.phi[1]) - model.Ue[1]*model.Um[1]*sin(model.phi[2]));
	}

	//amueboonecpv=+2.*ue(6)*um(6)*(-ue(4)*um(4)*sin(phi(2))-ue(5)*um(5)*sin(phi(3)))

	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));
	oscCon.aMuE[3] = -4.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*cos(model.phi[0]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[3] = 2.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*sin(model.phi[0]);
		else oscCon.aMuE_CPV[3] = -2.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*sin(model.phi[0]);
	}

	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));
	oscCon.aMuE[4] = -4.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*cos(model.phi[1]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[4] = 2.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*sin(model.phi[1]);
		else    oscCon.aMuE_CPV[4] = -2.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*sin(model.phi[1]);
	}

	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));
	oscCon.aMuE[5] = -4.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*cos(model.phi[2]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[5] = 2.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*sin(model.phi[2]);
		else    oscCon.aMuE_CPV[5] = -2.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*sin(model.phi[2]);
	}

	return oscCon;
}

oscContribution getOscContributionsNueDis(neutrinoModel model){
	oscContribution oscCon;

	oscCon.dm2[0] = pow(model.mNu[0],2);
	oscCon.aEE[0] = -4 * pow(model.Ue[0],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[1] = pow(model.mNu[1],2);
	oscCon.aEE[1] = -4 * pow(model.Ue[1],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[2] = pow(model.mNu[2],2);
	oscCon.aEE[2] = -4 * pow(model.Ue[2],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));
	oscCon.aEE[3] = -4. * pow(model.Ue[0],2) * pow(model.Ue[1],2);

	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));
	oscCon.aEE[4] = -4. * pow(model.Ue[0],2) * pow(model.Ue[2],2);

	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));
	oscCon.aEE[5] = -4. * pow(model.Ue[1],2) * pow(model.Ue[2],2);

	return oscCon;
}

oscContribution getOscContributionsNumuDis(neutrinoModel model){
	oscContribution oscCon;

	oscCon.aMuMu[0] = -4. * pow(model.Um[0],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[0] = pow(model.mNu[0],2);

	oscCon.aMuMu[1] = -4. * pow(model.Um[1],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[1] = pow(model.mNu[1],2);

	oscCon.aMuMu[2] = -4. * pow(model.Um[2],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[2] = pow(model.mNu[2],2);

	oscCon.aMuMu[3] = -4. * pow(model.Um[0],2) * pow(model.Um[1],2);
	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));

	oscCon.aMuMu[4] = -4. * pow(model.Um[0],2) * pow(model.Um[2],2);
	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));

	oscCon.aMuMu[5] = -4. * pow(model.Um[2],2) * pow(model.Um[1],2);
	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));

	return oscCon;
}

// Get your models sorted out (all of this has been tested to shit. we're good.)
neutrinoModel initializeMarkovParams(){

    // Initialize new model!
    neutrinoModel modelOld;
    modelOld.zero();
    bool reject;
    int noOfCPFactors = noOfSteriles*(noOfSteriles-1)/2;

    do{
        // Fill up our random array
        RanGen.RndmArray(13,ran);

        // Initial Params for Mass and Mixing
        if(noOfSteriles > 0){
            for(int i = 0; i < noOfSteriles; i++){
                modelOld.mNu[i] = pow(10., (TMath::Log10(dm2Min[i]) + ran[noOfSteriles*i]*TMath::Log10(dm2Max[i]/dm2Min[i]))/2);
                modelOld.Ue[i] = UMin + ran[3*i+1] * (UMax - UMin);
                modelOld.Um[i] = UMin + ran[3*i+2] * (UMax - UMin);
            }
        }
        // Now, let's do the CP factors, phi
        if(noOfCPFactors > 0){
            for(int i = 0; i < noOfCPFactors; i++){
                modelOld.phi[i] = double(ran[noOfSteriles*3+i]*2*TMath::Pi());
                if(CPConserving){
                    if(modelOld.phi[i] < TMath::Pi()) modelOld.phi[i] = 0;
                    else  modelOld.phi[i] = TMath::Pi();
                }
            }
        }
        reject = rejectModel(modelOld);
    }while(reject);

    return modelOld;
}


bool rejectModel(neutrinoModel model){

    int noOfCPFactors = noOfSteriles*(noOfSteriles-1)/2;
    // Now, we'll reject the model if matrix elements are too large
    reject1 = pow(model.Ue[0],2) + pow(model.Um[0],2) > UMaxSq || pow(model.Ue[1],2) + pow(model.Um[1],2) > UMaxSq || pow(model.Ue[2],2) + pow(model.Um[2],2) > UMaxSq || pow(model.Ue[0],2) + pow(model.Ue[1],2) + pow(model.Ue[2],2) > UMaxSq || pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2) > UMaxSq;

    // Another condition can be applied to avoid a negative under the square root for atmospheric neutrinos
    if(ATMOSPHERICProcess == 1){
        double dmuMax = .25;

        double A = (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2)) +
            pow((model.Um[0]*model.Um[1]),2) + pow((model.Um[0]*model.Um[2]),2) + pow((model.Um[1]*model.Um[2]),2);

        reject1 = reject1 || (1. - 4*A) < pow(1-2*dmuMax,2);
    }

    // More rejection conditions!
    reject2 = false;
    reject3 = false;
    reject4 = false;
    if(noOfSteriles > 1){
        // DEGENERACY FIX
        reject2 = model.mNu[1] < model.mNu[0] || abs(pow(model.mNu[1],2) - pow(model.mNu[0],2)) < dm2Min[2];
    }
    if(noOfSteriles > 2){
        // DEGENERACY FIX
        reject3 = model.mNu[2] < model.mNu[0] || model.mNu[2] < model.mNu[1] || abs(pow(model.mNu[2],2) - pow(model.mNu[0],2)) < dm2Min[2] || abs(pow(model.mNu[2],2) - pow(model.mNu[1],2)) < dm2Min[2];
    }

    if(scanType == 2){
        // For the Markov chain case, gotta check a few more things. Essentially whether or not we've stepped out of bounds.
        if(noOfSteriles > 0){
            for(int i = 0; i < noOfSteriles; i++){
                reject4 = reject4 || pow(model.mNu[i],2) < dm2Min[i] || pow(model.mNu[i],2) > dm2Max[i];

                if(usingUe)  reject4 = reject4 || model.Ue[i] < UMin || model.Ue[i] > UMax;
                if(usingUm)  reject4 = reject4 || model.Um[i] < UMin || model.Um[i] > UMax;
            }
        }
        if(noOfCPFactors > 0){
            for(int i = 0; i < noOfCPFactors; i++){
                reject4 = reject4 || model.phi[i] < 0 || model.phi[i] > 2*TMath::Pi();
            }
        }
    }

    return (reject1 || reject2 || reject3 || reject4);
}

bool jobOpt(){
    // Here, we're going to read out the jobOptions.txt file and assign those parameters.

    std::vector < double > paraVal;

    // Fill up paraVal vector
    std::string line;
    ifstream file;
    file.open(jobOptLoc+"jobOptions.txt");
     while (std::getline(file, line) ){
         std::istringstream is_line(line);
         std::string key;
         if( std::getline(is_line, key, '=')){
             std::string value;
             if(std::getline(is_line, value)){
                 paraVal.push_back(atof(value.c_str()));
             }
         }
    }

    // Assign those values.
    noOfSteriles =  paraVal[0];         UMax =          paraVal[1];     UMaxSq =            paraVal[2];
    CPConserving =  paraVal[3];         scanType =      paraVal[4];     gridPoints =        paraVal[5];
    jobID=          paraVal[6];         nMCGen=         paraVal[7];     rndInit=            paraVal[8];
    BugeyProcess=   paraVal[9];         CCFRProcess=    paraVal[10];    CDHSProcess=        paraVal[11];
    CHOOZProcess=   paraVal[12];        KARMENProcess=  paraVal[13];    LSNDProcess=        paraVal[14];
    NOMADProcess=   paraVal[15];        MBProcess=      paraVal[16];    MBProcessNubar=     paraVal[17];
    ATMOSPHERICProcess= paraVal[18];    NUMIProcess=    paraVal[19];    MINOSProcess=       paraVal[20];
    MINOSNCProcess= paraVal[21];        GALLIUMProcess= paraVal[22];    ReactorAnomaly=     paraVal[23];
    XSECProcess=    paraVal[24];        MBDISProcess=   paraVal[25];    MBDISProcessNubar=  paraVal[26];
    chi2Cut =       paraVal[27];        stepSize =      paraVal[28];    temperature =       paraVal[29];

    return true;
}

void bruteforce3p1(){

	globInit();
	globChisq(0);
    return;
}
#if !defined(__CINT__) || defined (__MAKECINT__)
int main()
{
	bruteforce3p1();
	return 0;
}
#endif
