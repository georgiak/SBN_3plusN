/* ------------------------------------------//
Created by Davio Cianci and Georgia Karagiorgi
Jan 25th, 2016

Notes:
 - found error in fortran code in ourpred_ws for minos

------------------------------------------// */

#include "globalFit.h"
#include <getopt.h>
#include "TStopwatch.h"

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

TRandom RanGen;

double chi2Log, chi2LogOld;

// Ntuple Variables
float chi2, dof, step, temp, gof, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56;

boonePackage mbNuPack, mbNubarPack; atmPackage atmPack; numiPackage numiPack; sinSqPackage lsndPack, karmenPack; galPackage galPack; cdhsPackage cdhsPack;
minosPackage minosPack; minosncPackage minosncPack; booneDisPackage mbNuDisPack, mbNubarDisPack; nomadPackage nomadPack; ccfrPackage ccfrPack;
bugeyPackage bugeyPack; choozPackage choozPack; xsecPackage xsecPack; booneDisPlusPackage mbNuDisPlusPack, mbNubarDisPlusPack; boonePlusPackage mbNuPlusPack, mbNubarPlusPack;

double ntupleFinish = 3.;
bool plusmode = false;

int globInit(){

    // read jobOption file and fill variables
    jobOpt();
    dm2Max[0] = 100.;   dm2Max[1] = 100.;    dm2Max[2] = 100.;
    dm2Min[0] = .01;     dm2Min[1] = .01;    dm2Min[2] = .01;
	UMax = .5;

    // INITIALIZATIONS
	std::cout << "Start initializations!" << std::endl;

    dm2VecInit(.01, 100.);
	if(MBProcess){
		if(!plusmode)	mbNuPack = mbNuInit();
		else	mbNuPlusPack = mbNuInitPlus();
	}
	if(MBProcessNubar){
		if(!plusmode)	mbNubarPack = mbNubarInit();
		else 	mbNubarPlusPack = mbNubarInitPlus();
	}
	if(ATMOSPHERICProcess) atmPack = atmInit();
	if(NUMIProcess) numiPack = numiInit();
	if(LSNDProcess) lsndPack = lsndInit();
	if(KARMENProcess) karmenPack = karmenInit();
	if(GALLIUMProcess) galPack = galInit();
	if(MINOSProcess) minosPack = minosInit();
	if(MINOSNCProcess) minosncPack = minosncInit();
	if(MBDISProcess){
		if(!plusmode)	mbNuDisPack = mbNuDisInit();
		else mbNuDisPlusPack = mbNuDisInitPlus();
	}
	if(MBDISProcessNubar){
		if(!plusmode)	mbNubarDisPack = mbNubarDisInit();
		else mbNubarDisPlusPack = mbNubarDisInitPlus();
	}
	if(NOMADProcess) nomadPack = nomadInit();
	if(CCFRProcess) ccfrPack = ccfrInit();
	if(CDHSProcess) cdhsPack = cdhsInit();
	if(BugeyProcess) bugeyPack = bugeyInit();
	if(CHOOZProcess) choozPack = choozInit();
	if(XSECProcess) xsecPack = xsecInit();

    getNDF();
	if(XSECProcess || BugeyProcess || CHOOZProcess)	myMinInit();

	if(debug) std::cout << "DOF: " << ndf << std::endl;

    std::cout << "Alright! Detector stuff successfully initialized!" << std::endl;

	return 1;
}

int globChisq(int ind){

	TStopwatch * watch = new TStopwatch();

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
	step = RanGen.Rndm() * stepSize;
	temp = RanGen.Rndm() * temperature;

    // Initialize the parameters we'll be using
    neutrinoModel nuModel;
    neutrinoModel nuModelOld;

    // Declare the chisq result
    chisqStruct chisqDetector, chisqTotal;

    std::cout << "Number of MC models = " << nMCGen << std::endl;

    // Initialize Markov Chain Parameters
    std::cout << "Initializing Markov chain parameters..." << std::endl;
    chi2Log = 0;    chi2LogOld = 0;
    nuModelOld = initializeMarkovParams();

    // CHI2 CALCULATIONS
    std::cout << "Start generating neutrino mass and mixing models..." << std::endl;
	int count = 0;
    for(iMCGen = 1; iMCGen <= nMCGen; iMCGen++){
        std::cout << "iMCGen = " << iMCGen << "/" << nMCGen << "... " << std::endl;

        if(count == 0){
			nuModel = initializeMarkovParams();
		}
		else
			nuModel = newModel(nuModelOld);

		if (count > nMCGen*ntupleFinish)
			break;
        chisqTotal.zero();  chisqDetector.zero();

		// For testing, we settin' our own parameters!
		//nuModel.zero();
		//nuModel.Ue[0] = .11; 	nuModel.Um[0] = .12;	nuModel.mNu[0] = sqrt(.9);
		//nuModel.Ue[1] = .11;	nuModel.Um[1] = .17; 	nuModel.mNu[1] = sqrt(17); 	nuModel.phi[0] = 1.6*TMath::Pi();
		//nuModel.Ue[2] = .11;	nuModel.Um[2] = .14;	nuModel.mNu[2] = sqrt(22);	nuModel.phi[1] = .28*TMath::Pi(); nuModel.phi[2] = 1.4*TMath::Pi();

		// Now, let's actually start calculating the chisq
		if(ATMOSPHERICProcess == 1){
			chisqDetector = getChi2Atm(nuModel, atmPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "ATM: " << chisqDetector.chi2 << std::endl;
		}
		if(MINOSProcess == 1){
            chisqDetector = getChi2Minos(nuModel, minosPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Minos: " << chisqDetector.chi2 << std::endl;
        }
		if(CCFRProcess == 1){
            chisqDetector = getChi2CCFR(nuModel, ccfrPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "CCFR: " << chisqDetector.chi2 << std::endl;
        }
        if(CDHSProcess == 1){
            chisqDetector = getChi2CDHS(nuModel, cdhsPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "CDHS: " << chisqDetector.chi2 << std::endl;
        }
		if(LSNDProcess == 1){
            chisqDetector = getLogLikelihood(nuModel, 5, lsndPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "LSND: " << chisqDetector.chi2 << std::endl;
        }
        if(KARMENProcess == 1){
            chisqDetector = getLogLikelihood(nuModel, 9, karmenPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Karmen: " << chisqDetector.chi2 << std::endl;
        }
        if(NUMIProcess == 1){
            chisqDetector = getChi2Numi(nuModel, numiPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Numi: " << chisqDetector.chi2 << std::endl;
        }
        if(NOMADProcess == 1){
            chisqDetector = getChi2Nomad(nuModel, nomadPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Nomad: " << chisqDetector.chi2 << std::endl;
        }
        if(MBProcess == 1){
			if(!plusmode)	chisqDetector = getChi2Boone(nuModel, mbNuPack, false);
			else chisqDetector = getChi2BoonePlus(nuModel, mbNuPlusPack, false);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
			chisqTotal.chi2 += chisqDetector.chi2;
            chisqTotal.chi2_det += chisqDetector.chi2_det;
			if(debug) std::cout << "MB: " << chisqDetector.chi2 << std::endl;
        }
		if(GALLIUMProcess == 1){
			chisqDetector = getChi2Gallium(nuModel, galPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Gal: " << chisqDetector.chi2 << std::endl;
		}
        if(MBProcessNubar == 1){
			if(!plusmode)	chisqDetector = getChi2Boone(nuModel, mbNubarPack, true);
			else chisqDetector = getChi2BoonePlus(nuModel, mbNubarPlusPack, true);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
            chisqTotal.chi2_det += chisqDetector.chi2_det;
			if(debug) std::cout << "MBNubar: " << chisqDetector.chi2 << std::endl;
        }
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(BugeyProcess == 1){
            chisqDetector = getChi2Bugey(nuModel, bugeyPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Bugey: " << chisqDetector.chi2 << std::endl;
        }
		if(chisqTotal.chi2 > chi2Cut) continue;
		if(XSECProcess == 1){
            chisqDetector = getChi2Xsec(nuModel, xsecPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "Xsec: " << chisqDetector.chi2 << std::endl;
        }
		if(chisqTotal.chi2 > chi2Cut) continue;

		if(MBDISProcessNubar == 1){
			if(!plusmode)	chisqDetector = getChi2MBDis(nuModel, mbNubarDisPack);
			else chisqDetector = getChi2MBDisPlus(nuModel, mbNubarDisPlusPack);
			if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
			chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "MBDisNubar: " << chisqDetector.chi2 << std::endl;
		}
		if(chisqTotal.chi2 > chi2Cut) continue;
        if(MBDISProcess == 1){
            if(!plusmode)	chisqDetector = getChi2MBDis(nuModel, mbNuDisPack);
			else chisqDetector = getChi2MBDisPlus(nuModel, mbNuDisPlusPack);
            if(chisqDetector.chi2 > chi2Cut || chisqDetector.chi2 < 0) continue;
            chisqTotal.chi2 += chisqDetector.chi2;
			if(debug) std::cout << "MBDis: " << chisqDetector.chi2 << std::endl;
        }

		if(chisqTotal.chi2 > chi2Cut) continue;

		chi2Log = abs(chisqTotal.chi2-trim);
		//gof = TMath::Prob(chi2Log,ndf);

		if(ran[12] < TMath::Min(1., exp(-(chi2Log - chi2LogOld)/temp)) || iMCGen == 1){
            chi2LogOld = chi2Log;
            nuModelOld = nuModel;
        }

        // Fill Ntuple
        chi2 = chisqTotal.chi2;
        dof = ndf;  m4 = nuModel.mNu[0];    m5 = nuModel.mNu[1];    m6 = nuModel.mNu[2];    ue4 = nuModel.Ue[0];    ue5 = nuModel.Ue[1];    ue6 = nuModel.Ue[2];
        um4 = nuModel.Um[0];    um5 = nuModel.Um[1];    um6 = nuModel.Um[2]; phi45 = nuModel.phi[0];    phi46 = nuModel.phi[1]; phi56 = nuModel.phi[2];
        chi2Nt->Fill(chi2, step, temp, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);

		count++;
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
            model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[noOfSteriles*i] - .5)*2*step*TMath::Log10(dm2Max[i]/dm2Min[i]))/2);

            if(usingUe)     model.Ue[i] = modelOld.Ue[i] + 2*(ran[noOfSteriles*i+1] - 0.5)*(UMax - UMin)*step;
            else    model.Ue[i] = 0.;

            if(usingUm)     model.Um[i] = modelOld.Um[i] + 2*(ran[noOfSteriles*i+2] - 0.5)*(UMax - UMin)*step;
            else    model.Um[i] = 0.;
        }
        if(noOfCPFactors > 0){
            for(int j = 0; j < noOfCPFactors; j++){
                model.phi[j] = modelOld.phi[j] + 2.*(ran[noOfSteriles*3 + j] - 0.5)*2.*TMath::Pi()*step;
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
    noOfSteriles =  	paraVal[0];		UMax =          	paraVal[1];		UMaxSq =        	paraVal[2];
	CPConserving =  	paraVal[3];		gridPoints =    	paraVal[4];		jobID=          	paraVal[5];
	nMCGen=         	paraVal[6];		rndInit=        	paraVal[7];	    BugeyProcess=   	paraVal[8];
	CCFRProcess=    	paraVal[9];		CDHSProcess=    	paraVal[10];	CHOOZProcess=   	paraVal[11];
	KARMENProcess=  	paraVal[12];	LSNDProcess=    	paraVal[13];    NOMADProcess=   	paraVal[14];
	MBProcess=      	paraVal[15];	MBProcessNubar= 	paraVal[16];    ATMOSPHERICProcess= paraVal[17];
	NUMIProcess=    	paraVal[18];	MINOSProcess=   	paraVal[19];	MINOSNCProcess= 	paraVal[20];
	GALLIUMProcess= 	paraVal[21];	ReactorAnomaly= 	paraVal[22];    XSECProcess=    	paraVal[23];
	MBDISProcess=   	paraVal[24];	MBDISProcessNubar=  paraVal[25];    chi2Cut =       	paraVal[26];
	stepSize =      	paraVal[27];	temperature =       paraVal[28];	trim =       		paraVal[29];
    return true;
}

#if !defined(__CINT__) || defined (__MAKECINT__)
int main(int argc, char* argv[])
{
	using namespace std;

	int index;
	int iarg = 0;
	opterr=1;
	const struct option longopts[] = {
	{"debug",	 		no_argument, 		0, 'D'},
	{"singlepoint",	 	no_argument, 		0, 's'},
	{"grid",	 		no_argument, 		0, 'g'},
	{"plus",			no_argument,		0, 'P'},
	};


	while(iarg != -1){
		iarg = getopt_long(argc,argv, "DsgP", longopts, &index);

		switch(iarg)
		{
			case 'D':
				debug = true;
				break;
			case 'g':
				jobOptLoc = "";
				dataLoc = "";
				break;
			case 'P':
				plusmode = true;
				break;
		}
	}

	globInit();
	globChisq(0);
	return 0;
}
#endif
