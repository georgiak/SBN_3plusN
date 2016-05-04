/* ------------------------------------------//
Created by Davio Cianci and Georgia Karagiorgi
Jan 25th, 2016

Notes:
 - found error in fortran code in ourpred_ws for minos

------------------------------------------// */

#include "globalFit.h"

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
float UMin = 0.01;

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

	jobOptLoc = ""; ///Users/dcianci/Physics/SBN_3plusN/GlobalFits/inputs/"; // /pnfs/lar1nd/scratch/users/dcianci/fits/";
	dataLoc = ""; ///Users/dcianci/Physics/SBN_3plusN/GlobalFits/data/"; ///pnfs/lar1nd/scratch/users/dcianci/fits/data";

    // read jobOption file and fill variables
    jobOpt();
    dm2Max[0] = 100.;   dm2Max[1] = 100.;    dm2Max[2] = 100.;
    dm2Min[0] = .1;     dm2Min[1] = .01;    dm2Min[2] = .01;

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

	// create the ntuple where the results are going to go
	std::string jid = Form("%d",(jobID + ind));
	//std::string outfile = target + "globFit_" + jid + ".root";
	std::string outfile = "globFit.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

	TNtuple *chi2Nt = new TNtuple("chi2Nt","chi2Nt","chi2:step:temp:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	//RanGen.SetSeed(rndInit + ind);
	RanGen.SetSeed(0);

	//We're gonna set a random step size and temp size for now - that way we can play with it and see what works best, yeah?
	step = RanGen.Rndm() * .8;
	temp = RanGen.Rndm() * 20.;

    // Initialize the parameters we'll be using
    neutrinoModel nuModel;
    neutrinoModel nuModelOld;

    // Initialize the chisq result
    chisqStruct chisqDetector, chisqTotal;

    std::cout << "Scantype = " << scanType << std::endl;
    if(scanType == 2)   std::cout << "Number of MC models = " << nMCGen << std::endl;
    if(scanType == 1)   std::cout << "Number of grid points = " << gridPoints << std::endl;

    // Initialize Markov Chain Parameters
    if(scanType == 2){
        std::cout << "Initializing Markov chain parameters..." << std::endl;

        chi2Log = 0;    chi2LogOld = 0;
        nuModelOld = initializeMarkovParams();
    }
	
    // CHI2 CALCULATIONS
    std::cout << "Start generating neutrino mass and mixing models..." << std::endl;

    // If we're doing a Markov Scan, here we are!
    if(scanType == 2){

        for(iMCGen = 1; iMCGen <= nMCGen; iMCGen++){
            std::cout << "iMCGen = " << iMCGen << "/" << nMCGen << "... " << std::endl;

            bool rejectCut; // calculated chisq fails cut
            do{
                rejectCut = true;
                if(iMCGen == 1)
					nuModel = initializeMarkovParams();
				else
					nuModel = newModel(nuModelOld);
                chisqTotal.zero();  chisqDetector.zero();

				// For testing, we settin' our own parameters!
				//nuModel.zero();
				//nuModel.Ue[0] = .11; 	nuModel.Um[0] = .12;	nuModel.mNu[0] = sqrt(.9);
				//nuModel.Ue[1] = .11;	nuModel.Um[1] = .17; 	nuModel.mNu[1] = sqrt(17); 	nuModel.phi[0] = 1.6*TMath::Pi();
				//nuModel.Ue[2] = .11;	nuModel.Um[2] = .14;	nuModel.mNu[2] = sqrt(22);	nuModel.phi[1] = .28*TMath::Pi(); nuModel.phi[2] = 1.4*TMath::Pi();

				//nuModel.zero();
				//nuModel.Ue[0] = .12; 	nuModel.Um[0] = .013;	nuModel.mNu[0] = sqrt(.92);
				//nuModel.Ue[1] = .16;	nuModel.Um[1] = .019; 	nuModel.mNu[1] = sqrt(7.2); 	nuModel.phi[0] = 1.6*TMath::Pi();
				//nuModel.Ue[2] = .069;	nuModel.Um[2] = .15;	nuModel.mNu[2] = sqrt(18);	nuModel.phi[1] = .28*TMath::Pi(); nuModel.phi[2] = 1.4*TMath::Pi();

                // Now, let's actually start calculating the chisq
                if(ATMOSPHERICProcess == 1){
                    chisqDetector = getChi2Atm(nuModel, atmPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "ATM: " << chisqDetector.chi2 << std::endl;
                }
                if(MBProcess == 1){
                    chisqDetector = getChi2Boone(nuModel, mbNuPack, false);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

					chisqTotal.chi2 += chisqDetector.chi2;
                    chisqTotal.chi2_det += chisqDetector.chi2_det;
					if(debug) std::cout << "MB: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(MBProcessNubar == 1){
                    chisqDetector = getChi2Boone(nuModel, mbNubarPack, true);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
                    chisqTotal.chi2_det += chisqDetector.chi2_det;
					if(debug) std::cout << "MBNubar: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(LSNDProcess == 1){
                    chisqDetector = getLogLikelihood(nuModel, 5, lsndPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "LSND: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(KARMENProcess == 1){
                    chisqDetector = getLogLikelihood(nuModel, 9, karmenPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Karmen: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(GALLIUMProcess == 1){
                    chisqDetector = getChi2Gallium(nuModel, galPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Gal: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(NUMIProcess == 1){
                    chisqDetector = getChi2Numi(nuModel, numiPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Numi: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(NOMADProcess == 1){
                    chisqDetector = getChi2Nomad(nuModel, nomadPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Nomad: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(MINOSProcess == 1){
                    chisqDetector = getChi2Minos(nuModel, minosPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Minos: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(MINOSNCProcess == 1){
                    // All of this has just been calculated elsewhere, so it's nice and quick.
                    double theta24_minosnc = asin(nuModel.Um[0]);
                    chisqDetector.chi2 = pow(theta24_minosnc - minosncPack.theta24_bestfit,2) / pow(minosncPack.sd_theta24,2);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Minosnc: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(MBDISProcess == 1){
                    chisqDetector = getChi2MBDis(nuModel, mbNuDisPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "MBDis: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(MBDISProcessNubar == 1){
                    chisqDetector = getChi2MBDis(nuModel, mbNubarDisPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "MBDisNubar: " << chisqDetector.chi2 << std::endl;
                }
                if(CCFRProcess == 1){
                    chisqDetector = getChi2CCFR(nuModel, ccfrPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "CCFR: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
                if(CDHSProcess == 1){
                    chisqDetector = getChi2CDHS(nuModel, cdhsPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "CDHS: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
				if(BugeyProcess == 1){
                    chisqDetector = getChi2Bugey(nuModel, bugeyPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Bugey: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
				if(CHOOZProcess == 1){
                    chisqDetector = getChi2Chooz(nuModel, choozPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Chooz: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
				if(XSECProcess == 1){
                    chisqDetector = getChi2Xsec(nuModel, xsecPack);
                    if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;

                    chisqTotal.chi2 += chisqDetector.chi2;
					if(debug) std::cout << "Xsec: " << chisqDetector.chi2 << std::endl;
                }
				if(chisqTotal.chi2 > chi2Cut) continue;
				rejectCut = false;
            }while(rejectCut);

            std::cout << "chisqtotal: " << chisqTotal.chi2 << std::endl;
			if(chisqTotal.chi2 < 0.) {
				std::cout << "WHOAWHOAWHOAAAH HOLD THE GODDAMN PHONE" << std::endl;
				return 0;
			}

            chi2Log = chisqTotal.chi2;
			gof = TMath::Prob(chi2Log,ndf);

			//for testing:
			temperature = temp;

            if(scanType == 2){
                if(ran[12] < TMath::Min(1., exp(-(chi2Log - chi2LogOld)/temperature)) || iMCGen == 1){
                    chi2LogOld = chi2Log;
                    nuModelOld = nuModel;
                }
            }

            // Fill Ntuple
            chi2 = chisqTotal.chi2;
            dof = ndf;  m4 = nuModel.mNu[0];    m5 = nuModel.mNu[1];    m6 = nuModel.mNu[2];    ue4 = nuModel.Ue[0];    ue5 = nuModel.Ue[1];    ue6 = nuModel.Ue[2];
            um4 = nuModel.Um[0];    um5 = nuModel.Um[1];    um6 = nuModel.Um[2]; phi45 = nuModel.phi[0];    phi46 = nuModel.phi[1]; phi56 = nuModel.phi[2];
            chi2Nt->Fill(chi2, step, temp, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);

			// For testing
			//break;

		}
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

neutrinoModel newModel(neutrinoModel modelOld){

	// for testing
	stepSize = step;

    // Initialize new model!
    neutrinoModel model;
    model.zero();
    bool reject;
    int noOfCPFactors = noOfSteriles*(noOfSteriles-1)/2;

    do{
        // Generate some random numbers!
        RanGen.RndmArray(13,ran);

        // Alright, let's step forward with these masses and mixing matrix elements!
        for(int i = 0; i < noOfSteriles; i++){
            model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[noOfSteriles*i] - .5)*2*stepSize*TMath::Log10(dm2Max[i]/dm2Min[i]))/2);

            if(usingUe)     model.Ue[i] = modelOld.Ue[i] + 2*(ran[noOfSteriles*i+1] - 0.5)*(UMax - UMin)*stepSize;
            else    model.Ue[i] = 0.;

            if(usingUm)     model.Um[i] = modelOld.Um[i] + 2*(ran[noOfSteriles*i+2] - 0.5)*(UMax - UMin)*stepSize;
            else    model.Um[i] = 0.;
        }
        if(noOfCPFactors > 0){
            for(int j = 0; j < noOfCPFactors; j++){
                model.phi[j] = modelOld.phi[j] + 2.*(ran[noOfSteriles*3 + j] - 0.5)*2.*TMath::Pi()*stepSize;
                if(CPConserving == 1){
                    if(model.phi[j] < TMath::Pi())    model.phi[j] = 0;
                    else model.phi[j] = TMath::Pi();
                }
            }
        }
        reject = rejectModel(model);
    }while(reject);

   return model;
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


void globalFit(){

	globInit();
	globChisq(0);
    return;
}
#if !defined(__CINT__) || defined (__MAKECINT__)
int main()
{
	globalFit();
	return 0;
}
#endif

// This is just a temp thing to make a plot for a talk
void makeMBPlot(){

	double binso[12] = {200.,  300.,  375.,  475.,  550.,  675.,  800.,  950.,  1100.,  1300.,  1500.,  3000.};
	double bins = 11;

	TH1D* hist_obs = new TH1D("hist_obs", "", bins,binso);
	TH1D* hist_pred_31 = new TH1D("hist_pred_31", "", bins,binso);
	TH1D* hist_pred_32 = new TH1D("hist_pred_32", "", bins,binso);
	TH1D* hist_pred_33 = new TH1D("hist_pred_33", "", bins,binso);

	TCanvas* c9 = new TCanvas("c3","Signal vs Prediction",700,700);
	c9->SetLeftMargin(.15);
	c9->SetBottomMargin(.15);
	c9->SetTopMargin(.05);
	c9->SetRightMargin(.05);
	c9->cd();

	hist_obs->SetTitle("Observed vs. Predicted Signal for MiniBooNE;Energy (MeV);Events/Bin");
	hist_obs->GetXaxis()->SetTitleOffset(1);
	hist_obs->GetYaxis()->SetTitleOffset(1.4);
	hist_obs->GetXaxis()->SetTitleFont(62);
	hist_obs->GetYaxis()->SetTitleFont(62);
	hist_obs->GetYaxis()->CenterTitle();
	hist_obs->GetXaxis()->CenterTitle();
	hist_obs->GetXaxis()->SetTitleSize(0.05);
	//hist_obs->GetXaxis()->SetLabelSize(0.04);
	hist_obs->GetXaxis()->SetLabelOffset(0.001);
	hist_obs->GetYaxis()->SetTitleSize(0.05);
	//hist_obs->GetYaxis()->SetLabelSize(0.04);
	hist_obs->SetStats(kFALSE);

	globInit();
	boonePackage pack = mbNuPack;

	int nBins = 11;
	int nBins_mu = 8;
	int nFOscEvts = 17204;

	neutrinoModel threep1, threep2, threep3;
	threep1.zero();
	threep1.Ue[0] = .15; 	threep1.Um[0] = .17;	threep1.mNu[0] = sqrt(.92);
	threep2.zero();
	threep2.Ue[0] = .15; 	threep2.Um[0] = .13;	threep2.mNu[0] = sqrt(.92);
	threep2.Ue[1] = .069;	threep2.Um[1] = .16; 	threep2.mNu[1] = sqrt(17.); 	threep2.phi[0] = 1.8*TMath::Pi();
	threep3.zero();
	threep3.Ue[0] = .11; 	threep3.Um[0] = .12;	threep3.mNu[0] = sqrt(.9);
	threep3.Ue[1] = .11;	threep3.Um[1] = .17; 	threep3.mNu[1] = sqrt(17.); 	threep3.phi[0] = 1.6*TMath::Pi();
	threep3.Ue[2] = .11;	threep3.Um[2] = .14;	threep3.mNu[2] = sqrt(22.);	threep3.phi[1] = .28*TMath::Pi(); threep3.phi[2] = 1.4*TMath::Pi();

	// Initialize contributions from the oscillation probability
    oscContribution oscCont;

	_signal.resize(nBins);
	_fullData.resize(nBins + nBins_mu);
	_prediction.resize(nBins + nBins_mu);

	double minEBins[nBins];
	double maxEBins[nBins];
	double ETru[nFOscEvts], LTru[nFOscEvts];

	// THREEP1
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

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep1, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_31->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
	    }

	// threep2// Initialize vars
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

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep2, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_32->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
	    }

	// threep3// Initialize vars
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

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep3, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_33->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
			hist_obs->SetBinContent(iB+1,_fullData[iB]/(binso[iB+1]-binso[iB]));
			hist_obs->SetBinError(iB+1,sqrt(_fullData[iB])/(binso[iB+1]-binso[iB]));
	    }

		// Okay, now draw it.
		hist_obs->SetMarkerStyle(20);
		hist_obs->SetMarkerColor(kBlack);
		hist_pred_33->SetLineWidth(3);
		hist_pred_33->SetLineColor(kMagenta);
		hist_pred_32->SetLineWidth(3);
		hist_pred_32->SetLineColor(kBlue);
		hist_pred_31->SetLineWidth(3);
		hist_pred_31->SetLineColor(kCyan);

		TLegend* legS=new TLegend(0.7,0.7,0.95,0.9);
		legS->SetFillStyle(0);
		legS->SetFillColor(0);
		legS->SetBorderSize(0);
		legS->SetTextFont(62);
		legS->SetTextSize(0.03);
		legS->AddEntry(hist_obs,"Observed","pe");
		legS->AddEntry(hist_pred_31,"3+1 Pred","f");
		legS->AddEntry(hist_pred_32,"3+2 Pred","f");
		legS->AddEntry(hist_pred_33,"3+3 Pred","f");

		hist_obs->Draw("p");
		hist_pred_31->Draw("h same");
		hist_pred_32->Draw("h same");
		hist_pred_33->Draw("h same");

		legS->Draw();

		c9->SaveAs("obsSignal.png");

}
