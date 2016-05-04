/* ------------------------------------------//
Created by Davio Cianci
Jan 25th, 2016

This takes all the different markov chains and smashes them together, finds the minimum chi2 and
then outputs two final ntuples that contain all points within the 90 and 99% CL

------------------------------------------// */

#include "TLegend.h"
#include "globalFit.h"
#include "TCut.h"
bool procOpt();

float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
float m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;
int steriles, nRuns;
std::string dataset, location, output;
std::string procOptLoc;

int ntProcess(){

	procOptLoc = "/lar1nd/app/users/dcianci/SBN_3plusN/GlobalFits/inputs/";
	procOpt();

	// Create tchain by linking all root ntuples
    std::cout << "Loading ntuple files..." << std::endl;
	float chi2min = 3000.;
    TChain * in_chain = new TChain("chi2Nt");

	for(int i = 0; i < nRuns; i++){
        std::string jid = Form("%i",i);
        std::string infile = "/pnfs/lar1nd/scratch/users/dcianci/" + location + "/globFit_" + jid + ".root";
        //std::cout << "Output File 1: " << infile << std::endl;
        TString inputFile = infile;
        in_chain->Add(inputFile);
    }

    in_chain->SetBranchAddress("chi2",&chi2);
	in_chain->SetBranchAddress("m4",&m4);
	in_chain->SetBranchAddress("ue4",&ue4);
	in_chain->SetBranchAddress("um4",&um4);
	in_chain->SetBranchAddress("m5",&m5);
	in_chain->SetBranchAddress("ue5",&ue5);
	in_chain->SetBranchAddress("um5",&um5);
	in_chain->SetBranchAddress("m6",&m6);
	in_chain->SetBranchAddress("ue6",&ue6);
	in_chain->SetBranchAddress("um6",&um6);
	in_chain->SetBranchAddress("phi45",&phi45);
	in_chain->SetBranchAddress("phi46",&phi46);
	in_chain->SetBranchAddress("phi56",&phi56);

    std::cout << in_chain->GetEntries() << std::endl;

    // Find chi2Min
    for(int i = 0; i < in_chain->GetEntries(); i++){
        in_chain->GetEntry(i);
		if(chi2 < chi2min){
			chi2min = chi2;
        	m4_min = m4;
			ue4_min = ue4;
			um4_min = um4;
			m5_min = m5;
			ue5_min = ue5;
			um5_min = um5;
			m6_min = m6;
			ue6_min = ue6;
			um6_min = um6;
			phi45_min = phi45;
			phi46_min = phi46;
			phi56_min = phi56;
		}
    }
    std::cout << chi2min << std::endl;
	std::cout << "m4_min: " << m4_min << std::endl;
	std::cout << "ue4_min: " << ue4_min << std::endl;
	std::cout << "um4_min: " << um4_min << std::endl;
	std::cout << "m5_min: " << m5_min << std::endl;
	std::cout << "ue5_min: " << ue5_min << std::endl;
	std::cout << "um5_min: " << um5_min << std::endl;
	std::cout << "m6_min: " << m6_min << std::endl;
	std::cout << "ue6_min: " << ue6_min << std::endl;
	std::cout << "um6_min: " << um6_min << std::endl;
	std::cout << "phi45_min: " << phi45_min << std::endl;
	std::cout << "phi46_min: " << phi46_min << std::endl;
	std::cout << "phi56_min: " << phi56_min << std::endl;

	// Make new file to fill with confidence level ntuples
	std::string jid = Form("/nt_3%i_",steriles);
	std::string outfile = output + jid + dataset + ".root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}
	TNtuple *chi2_99 = new TNtuple("chi2_99","chi2_99","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	TNtuple *chi2_90 = new TNtuple("chi2_90","chi2_90","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	for(int i = 0; i < in_chain->GetEntries(); i++){
        in_chain->GetEntry(i);

		if(chi2 - chi2min < 9.201)
			chi2_99->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
		if(chi2 - chi2min < 4.605)
			chi2_90->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
    }

	std::cout << "90CL has " << chi2_90->GetEntries() << " entries" << std::endl;
	std::cout << "99CL has " << chi2_99->GetEntries() << " entries" << std::endl;

	// Save Ntuple to file
	chi2_99->Write();
	chi2_90->Write();
	f->Close();

	// 90% 		4.605
	// 99%		9.201

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
    ntProcess();
    return 0;
}
# endif

void ntupleProcess(){
    ntProcess();
    return;
}
