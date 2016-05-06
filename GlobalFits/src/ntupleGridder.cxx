/* ------------------------------------------//
Created by Davio Cianci
Jan 25th, 2016

This takes all the different markov chains and smashes them together, finds the minimum chi2 and
then outputs two final ntuples that contain all points within the 90 and 99% CL

------------------------------------------// */

#include "TLegend.h"
#include "globalFit.h"
#include "TH3D.h"
bool procOpt();

float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
float m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;
int steriles, nRuns, type, raster;
std::string dataset, location, output;
std::string procOptLoc;

TNtuple * chi2_99, * chi2_90, * chi2_95;
int gridPoints;
float dmmin, dmmax, umin, umax;

int ntGridder(){

	procOptLoc = "/lar1nd/app/users/dcianci/SBN_3plusN/GlobalFits/inputs/";
	procOpt();

	gridPoints = 100;
	dmmin = 0.1;	dmmax = 10.;
	umin = 0;		umax = .5;

	// Grab processed ntuple
    std::cout << "Loading ntuple file..." << std::endl;
	std::string jid = Form("/nt_3%i_",steriles);
	std::string infile = output + jid + dataset + ".root";
	std::cout << "Input File: " << infile << std::endl;
	TString inputFile = infile;
	TFile *inf = new TFile(inputFile);
	TNtuple *chi2_99 = (TNtuple*)inf->Get("chi2_99");
	TNtuple *chi2_90 = (TNtuple*)inf->Get("chi2_90");


	// Make new file to fill with processed ntuples
	jid = Form("/nt_3%i_",steriles);
	std::string outfile = output + jid + dataset + "_processed.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

	TNtuple * chi2_99_pr = new TNtuple("chi2_99_pr","chi2_99_pr","m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	TNtuple * chi2_90_pr = new TNtuple("chi2_90_pr","chi2_90_pr","m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	if(steriles == 1){
		// Okay, this one's 'easy' using built-in root data structures.
		// Fill up a histogram so everything's in proper order
		std::cout << "Filling rastergram histo" << std::endl;
		TH3D * rastergram_99 = new TH3D("rg","rg",gridPoints,.1,10.,gridPoints,0,.5,gridPoints,0,.5);
		TH3D * rastergram_90 = new TH3D("rg","rg",gridPoints,.1,10.,gridPoints,0,.5,gridPoints,0,.5);

		// We'll deal with 99% first, yeah?
		chi2_99->SetBranchAddress("chi2",&chi2);
		chi2_99->SetBranchAddress("m4",&m4);
		chi2_99->SetBranchAddress("ue4",&ue4);
		chi2_99->SetBranchAddress("um4",&um4);

		for(int i = 0; i < chi2_99->GetEntries(); i++){
        	chi2_99->GetEntry(i);
			float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
			float ustep = (umax-umin)/float(gridPoints);
			float _m4 = ceil(TMath::Log10(m4/dmmin)/mstep);
			float _ue4 = ceil(ue4/ustep);
			float _um4 = ceil(um4/ustep);

			rastergram_99->SetBinContent(_m4,_ue4,_um4,chi2);
		}

		chi2_90->SetBranchAddress("chi2",&chi2);
		chi2_90->SetBranchAddress("m4",&m4);
		chi2_90->SetBranchAddress("ue4",&ue4);
		chi2_90->SetBranchAddress("um4",&um4);

		for(int i = 0; i < chi2_90->GetEntries(); i++){
        	chi2_90->GetEntry(i);
			float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
			float ustep = (umax-umin)/float(gridPoints);
			float _m4 = ceil(TMath::Log10(m4/dmmin)/mstep);
			float _ue4 = ceil(ue4/ustep);
			float _um4 = ceil(um4/ustep);

			rastergram_90->SetBinContent(_m4,_ue4,_um4,chi2);
		}

		// Now, fill the new ntuples
		for(int _m4 = 1; _m4 <= gridPoints; _m4++){
			for(int _ue4 = 1; _ue4 <= gridPoints; _ue4++){
				for(int _um4 = 1; _um4 <= gridPoints; _um4++){
					if(rastergram_99->GetBinContent(_m4,_ue4,_um4) > 0)
						chi2_99_pr->Fill(_m4,_ue4,_um4,0,0,0,0,0,0,0,0,0);
					if(rastergram_90->GetBinContent(_m4,_ue4,_um4) > 0)
						chi2_90_pr->Fill(_m4,_ue4,_um4,0,0,0,0,0,0,0,0,0);
				}
			}
		}
	}

	// Now, let's see how much this is really doing.
	std::cout << "For 99%%, we reduced " << chi2_99->GetEntries() << " events to " << chi2_99_pr->GetEntries() << std::endl;
	std::cout << "For 90%%, we reduced " << chi2_90->GetEntries() << " events to " << chi2_90_pr->GetEntries() << std::endl;

	// Save Ntuple to file
	chi2_99_pr->Write();
	chi2_90_pr->Write();
	f->Close();

    return 0;
}

bool procOpt(){
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
	std::istringstream is_line2(line);
	std::getline(is_line2,key,'=');
	std::getline(is_line2,value);
	nRuns = atoi(value.c_str());

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
