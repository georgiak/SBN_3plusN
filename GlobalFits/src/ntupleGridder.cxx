/* ------------------------------------------//
Created by Davio Cianci
Jan 25th, 2016

This takes in the processed 90% and 99% CL ntuples and processes them further, ridding them of redundancy
to make things easier for Mark. We split our space into a grid and check to see if we have multiple points that
more-or-less overlap.

------------------------------------------// */
#include "globalFit.h"
std::vector <int> getGridPoint(float m4, float ue4, float um4, float m5, float ue5, float um5, float m6, float ue6, float um6, float phi45, float phi46, float phi56);
bool sameGridPoint(std::vector <int> myGrid, float m4, float ue4, float um4, float m5, float ue5, float um5, float m6, float ue6, float um6, float phi45, float phi46, float phi56);
float getPtLin(int y, float min, float max);
float getPtLog(int y, float min, float max);
bool procOpt();

float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
int steriles, nRuns;
std::string dataset, location, output;
std::string procOptLoc;

double massMin = .01;
double massMax = 100.;
double UMin = 0.;
double UMax = .5;
double sigmaTol = 5.;		// For multiple points in the same grid, the max sigma we'll allow before we decide we must go finer
int nGridPoints = 100;

int ntGridder(){

	procOptLoc = "../inputs/";
	procOpt();

	// Make new file to fill with processed ntuples
	std::string jid = Form("/nt_3%i_",steriles);
	std::string outfile = output + jid + dataset + "_processed.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}
	TNtuple *chi2_99_pr = new TNtuple("chi2_99_pr","chi2_99_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	TNtuple *chi2_90_pr = new TNtuple("chi2_90_pr","chi2_90_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	// This method is either stupid or genius
	std::vector <int> entry90;
	std::vector <int> entry99;
	std::vector <int> mnmGrid;	// vector with mass-and-mixing parameter grid points.

	// Load in the ntuples
	std::string infile = output + jid + dataset + ".root";
	TString inputFile = infile;
	TFile *inf = new TFile(inputFile);
	TNtuple *chi_99 = (TNtuple*)inf->Get("chi2_99");
	TNtuple *chi_90 = (TNtuple*)inf->Get("chi2_90");

	// Initialize dumb index vector;
	for(int i = 0; i < chi_90->GetEntries(); i++)
		entry90.push_back(i);
	for(int i = 0; i < chi_99->GetEntries(); i++)
		entry99.push_back(i);

	// We have to do one at a time, so be patient.
	std::cout << "Processing 99CL ntuple..." << std::endl;
	std::cout << "It's got " << chi_99->GetEntries() << " entries, so be patient." << std::endl;

	// Assign branches
	chi_99->SetBranchAddress("chi2",&chi2);
	chi_99->SetBranchAddress("m4",&m4);
	chi_99->SetBranchAddress("ue4",&ue4);
	chi_99->SetBranchAddress("um4",&um4);
	chi_99->SetBranchAddress("m5",&m5);
	chi_99->SetBranchAddress("ue5",&ue5);
	chi_99->SetBranchAddress("um5",&um5);
	chi_99->SetBranchAddress("m6",&m6);
	chi_99->SetBranchAddress("ue6",&ue6);
	chi_99->SetBranchAddress("um6",&um6);
	chi_99->SetBranchAddress("phi45",&phi45);
	chi_99->SetBranchAddress("phi46",&phi46);
	chi_99->SetBranchAddress("phi56",&phi56);

	do{
		chi_99->GetEntry(entry99[0]);
		float chi2Mean = 0.;

		// Find it's n-dimensional grid point in parameter space
		mnmGrid = getGridPoint(m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);

		// Find all candidates within this grid point;
		std::vector <int> candyInd;		candyInd.resize(0);
		std::vector <int> candyChi2;	candyChi2.resize(0);
		candyInd.push_back(0);	candyChi2.push_back(chi2);

		for(int k = 1; k < entry99.size(); k++){
			chi_99->GetEntry(k);
			if(sameGridPoint(mnmGrid, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56)){
				candyInd.push_back(k);
				candyChi2.push_back(chi2);
			}
		}

		// Sweet, now find the mean and sigma chi2
		for(int i = 0; i < candyChi2.size(); i++){
			chi2Mean += candyChi2[i]/candyChi2.size();
		}
		float sqDev = 0.;
		for(int i = 0; i < candyChi2.size(); i++){
			sqDev += pow(candyChi2[i] - chi2Mean,2)/candyChi2.size();
		}
		if(sqrt(sqDev) > sigmaTol){
			std::cout << "Sigma's too big! It's : " << sqrt(sqDev) << std::endl;
			return 1;
		}

		// Now fill up the new ntuple
		chi2_99_pr->Fill(chi2Mean,getPtLog(mnmGrid[0],massMin,massMax),getPtLin(mnmGrid[1],UMin,UMax),
					getPtLin(mnmGrid[2],UMin,UMax),getPtLog(mnmGrid[3],massMin,massMax),getPtLin(mnmGrid[4],UMin,UMax),
					getPtLin(mnmGrid[5],UMin,UMax),getPtLog(mnmGrid[6],massMin,massMax),getPtLin(mnmGrid[7],UMin,UMax),
					getPtLin(mnmGrid[8],UMin,UMax),getPtLin(mnmGrid[9],0,TMath::Pi()),getPtLin(mnmGrid[10],0,TMath::Pi()),
					getPtLin(mnmGrid[11],0,TMath::Pi()));

		// Kill off entries we've used
		for(int i = 0; i < candyInd.size(); i++){
			entry99.erase(entry99.begin() + candyInd[i]);
		}

	}while(entry99.size() > 0);

	std::cout << "99CL ntuple processed! Now we have: " << chi2_99_pr->GetEntries() << " entries!" << std::endl;

	// - - - - - - - - - - - - - -

	// Now, the 90% CL
	std::cout << "Processing 90CL ntuple..." << std::endl;
	std::cout << "It's got " << chi_90->GetEntries() << " entries, so be patient." << std::endl;

	// Assign branches
	chi_90->SetBranchAddress("chi2",&chi2);
	chi_90->SetBranchAddress("m4",&m4);
	chi_90->SetBranchAddress("ue4",&ue4);
	chi_90->SetBranchAddress("um4",&um4);
	chi_90->SetBranchAddress("m5",&m5);
	chi_90->SetBranchAddress("ue5",&ue5);
	chi_90->SetBranchAddress("um5",&um5);
	chi_90->SetBranchAddress("m6",&m6);
	chi_90->SetBranchAddress("ue6",&ue6);
	chi_90->SetBranchAddress("um6",&um6);
	chi_90->SetBranchAddress("phi45",&phi45);
	chi_90->SetBranchAddress("phi46",&phi46);
	chi_90->SetBranchAddress("phi56",&phi56);

	do{
		chi_90->GetEntry(entry90[0]);
		float chi2Mean = 0.;

		// Find it's n-dimensional grid point in parameter space
		mnmGrid = getGridPoint(m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);

		// Find all candidates within this grid point;
		std::vector <int> candyInd;		candyInd.resize(0);
		std::vector <int> candyChi2;	candyChi2.resize(0);
		candyInd.push_back(0);	candyChi2.push_back(chi2);

		for(int k = 1; k < entry90.size(); k++){
			chi_90->GetEntry(k);
			if(sameGridPoint(mnmGrid, m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56)){
				candyInd.push_back(k);
				candyChi2.push_back(chi2);
			}
		}

		// Sweet, now find the mean and sigma chi2
		for(int i = 0; i < candyChi2.size(); i++){
			chi2Mean += candyChi2[i]/candyChi2.size();
		}
		float sqDev = 0.;
		for(int i = 0; i < candyChi2.size(); i++){
			sqDev += pow(candyChi2[i] - chi2Mean,2)/candyChi2.size();
		}
		if(sqrt(sqDev) > sigmaTol){
			std::cout << "Sigma's too big! It's : " << sqrt(sqDev) << std::endl;
			return 1;
		}

		// Now fill up the new ntuple
		chi2_90_pr->Fill(chi2,getPtLog(mnmGrid[0],massMin,massMax),getPtLin(mnmGrid[1],UMin,UMax),
					getPtLin(mnmGrid[2],UMin,UMax),getPtLog(mnmGrid[3],massMin,massMax),getPtLin(mnmGrid[4],UMin,UMax),
					getPtLin(mnmGrid[5],UMin,UMax),getPtLog(mnmGrid[6],massMin,massMax),getPtLin(mnmGrid[7],UMin,UMax),
					getPtLin(mnmGrid[8],UMin,UMax),getPtLin(mnmGrid[9],0,TMath::Pi()),getPtLin(mnmGrid[10],0,TMath::Pi()),
					getPtLin(mnmGrid[11],0,TMath::Pi()));

		// Kill off entries we've used
		for(int i = 0; i < candyInd.size(); i++){
			entry90.erase(entry90.begin() + candyInd[i]);
		}

	}while(entry90.size() > 0);

	std::cout << "90CL ntuple processed! Now we have: " << chi2_90_pr->GetEntries() << " entries!" << std::endl;

	chi2_90_pr->Write();
	chi2_99_pr->Write();
	f->Close();

	return 0;
}

std::vector <int> getGridPoint(float m4, float ue4, float um4, float m5, float ue5, float um5, float m6, float ue6, float um6, float phi45, float phi46, float phi56){
	// This function returns a vector with the coordinates in 12 dimensional space where our point lives.
	std::vector <int> gridpt;

	float massStep = TMath::Log10(massMax/massMin)/(nGridPoints+1);
	float mixingStep = (UMax - UMin)/(nGridPoints+1);
	float cpvStep = 2*TMath::Pi()/(nGridPoints+1);
	gridpt.push_back(floor(m4/massStep) - TMath::Log10(massMin));
	gridpt.push_back(floor(ue4/mixingStep));
	gridpt.push_back(floor(um4/mixingStep));
	gridpt.push_back(floor(m5/massStep) - TMath::Log10(massMin));
	gridpt.push_back(floor(ue5/mixingStep));
	gridpt.push_back(floor(um5/mixingStep));
	gridpt.push_back(floor(m6/massStep) - TMath::Log10(massMin));
	gridpt.push_back(floor(ue6/mixingStep));
	gridpt.push_back(floor(um6/mixingStep));
	gridpt.push_back(floor(phi45/cpvStep));
	gridpt.push_back(floor(phi46/cpvStep));
	gridpt.push_back(floor(phi56/cpvStep));

	return gridpt;
}

bool sameGridPoint(std::vector <int> myGrid, float m4, float ue4, float um4, float m5, float ue5, float um5, float m6, float ue6, float um6, float phi45, float phi46, float phi56){

	std::vector <int> gridpt;
	gridpt = getGridPoint(m4, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45, phi46, phi56);

	if(myGrid[0] != m4) return false;
	if(myGrid[1] != ue4) return false;
	if(myGrid[2] != um4) return false;
	if(myGrid[3] != m5) return false;
	if(myGrid[4] != ue5) return false;
	if(myGrid[5] != um5) return false;
	if(myGrid[6] != m6) return false;
	if(myGrid[7] != ue6) return false;
	if(myGrid[8] != um6) return false;
	if(myGrid[9] != phi45) return false;
	if(myGrid[10] != phi46) return false;
	if(myGrid[11] != phi56) return false;

	return true;
}

float getPtLin(int y, float min, float max){
	return y * (max - min)/(nGridPoints+1);
}
float getPtLog(int y, float min, float max){
	return pow((y + TMath::Log10(min)) * (TMath::Log10(max/min)/(nGridPoints+1)),20);
}

    // Here, we're going to read out the procOpt.txt file and assign those parameters.


    // Fill up paraVal vector
    std::string line;
    ifstream file;
    file.open(procOptLoc+"processOptions.txt");

	std::string key;
	std::string value;

	std::getline(file,line);
	std::istringstream is_line1(line);
	std::getline(is_line1,key,'=');
	std::getline(is_line1,value);
	steriles = atoi(value.c_str());

	std::getline(file,line);

	std::getline(file,line);
	std::istringstream is_line3(line);
	std::getline(is_line3,key,'=');
	std::getline(is_line3,dataset);

	std::getline(file,line);
	std::istringstream is_line4(line);
	std::getline(is_line4,key,'=');
	std::getline(is_line4,location);

	std::getline(file,line);
	std::istringstream is_line5(line);
	std::getline(is_line5,key,'=');
	std::getline(is_line5,output);

    return true;
}

#ifndef __CINT__
int main()
{
    ntGridder();
    return 0;
}
# endif

void ntupleGridder(){
    ntGridder();
    return;
}
