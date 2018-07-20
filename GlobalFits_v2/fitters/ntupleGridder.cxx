/* ------------------------------------------//
Created by Davio Cianci
Jan 25th, 2016

This takes all the different markov chains and smashes them together, finds the minimum chi2 and
then outputs two final ntuples that contain all points within the 90 and 99% CL

------------------------------------------// */

#include "TLegend.h"
#include "globalFit.h"
#include "TH3D.h"
#include <stdio.h>
bool procOpt();

std::string plotOutput = "plots";
float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
float chi2min, m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;
int steriles, nRuns, type, raster;
std::string dataset, location, output;
std::string procOptLoc;

TNtuple * chi2_99, * chi2_90, * chi2_95;
TNtuple * chi2_99_pr, * chi2_90_pr;
int gridPoints;
float dmmin, dmmax, umin, umax;

bool stringCompare(const string & l, const string & r)
{
	// We only want to compare the points, not the chisq
	return (l.substr(0,24) == r.substr(0,24));
}

int ntGridder(){

	procOptLoc = "inputs/";
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

	chi2_99->SetBranchAddress("chi2",&chi2);
	chi2_99->SetBranchAddress("m4",&m4);
	chi2_99->SetBranchAddress("ue4",&ue4);
	chi2_99->SetBranchAddress("um4",&um4);
	chi2_99->SetBranchAddress("m5",&m5);
	chi2_99->SetBranchAddress("ue5",&ue5);
	chi2_99->SetBranchAddress("um5",&um5);
	chi2_99->SetBranchAddress("m6",&m6);
	chi2_99->SetBranchAddress("ue6",&ue6);
	chi2_99->SetBranchAddress("um6",&um6);
	chi2_99->SetBranchAddress("phi45",&phi45);
	chi2_99->SetBranchAddress("phi46",&phi46);
	chi2_99->SetBranchAddress("phi56",&phi56);

    std::cout << chi2_99->GetEntries() << std::endl;

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

	double cl90, cl99;
	if(steriles == 1){
		cl90 = 6.25;
		cl99 = 11.34;
	}
	if(steriles == 2){
		cl90 = 12.02;
		cl99 = 18.48;
	}
	if(steriles == 3){
		cl90 = 18.55;
		cl99 = 26.22;
	}

	chi2_99_pr = new TNtuple("chi2_99_pr","chi2_99_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	chi2_90_pr = new TNtuple("chi2_90_pr","chi2_90_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	// Now, we're going to convert to a string of coordinates.
	std::cout << "Convert entries to strings of coords" << std::endl;

	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	float ustep = (umax-umin)/float(gridPoints);
	float phistep = 2*TMath::Pi()/float(gridPoints);
	std::vector<std::string> coords;
	coords.resize(chi2_99->GetEntries());

	int _m4;	int _m5;	int _m6;	int _ue4;	int _ue5;	int _ue6;	int _um4;	int _um5;	int _um6;	int _phi45;
	int _phi46;	int _phi56;
	std::string _m4s, _m5s, _m6s, _ue4s, _ue5s, _ue6s, _um4s, _um5s, _um6s, _phi45s, _phi46s, _phi56s, _chi2s;
	for(int i = 0; i < chi2_99->GetEntries(); i++){
		chi2_99->GetEntry(i);

		// Now, the chi2!
		stringstream stream;
		stream << fixed << setprecision(2) << chi2;
		_chi2s = stream.str();

		_m4 = TMath::Max(int(floor(TMath::Log10(m4/dmmin)/mstep)),0);
		if(_m4 < 10) _m4s = "0" + to_string(_m4);	else _m4s = to_string(_m4);
		_m5 = TMath::Max(int(floor(TMath::Log10(m5/dmmin)/mstep)),0);
		if(_m5 < 10) _m5s = "0" + to_string(_m5);	else _m5s = to_string(_m5);
		_m6 = TMath::Max(int(floor(TMath::Log10(m6/dmmin)/mstep)),0);
		if(_m6 < 10) _m6s = "0" + to_string(_m6);	else _m6s = to_string(_m6);
		_ue4 = TMath::Max(int(floor(ue4/ustep)),0);
		if(_ue4 < 10) _ue4s = "0" + to_string(_ue4);	else _ue4s = to_string(_ue4);
		_ue5 = TMath::Max(int(floor(ue5/ustep)),0);
		if(_ue5 < 10) _ue5s = "0" + to_string(_ue5);	else _ue5s = to_string(_ue5);
		_ue6 = TMath::Max(int(floor(ue6/ustep)),0);
		if(_ue6 < 10) _ue6s = "0" + to_string(_ue6);	else _ue6s = to_string(_ue6);
		_um4 = TMath::Max(int(floor(um4/ustep)),0);
		if(_um4 < 10) _um4s = "0" + to_string(_um4);	else _um4s = to_string(_um4);
		_um5 = TMath::Max(int(floor(um5/ustep)),0);
		if(_um5 < 10) _um5s = "0" + to_string(_um5);	else _um5s = to_string(_um5);
		_um6 = TMath::Max(int(floor(um6/ustep)),0);
		if(_um6 < 10) _um6s = "0" + to_string(_um6);	else _um6s = to_string(_um6);
		_phi45 = TMath::Max(int(floor(phi45/phistep)),0);
		if(_phi45 < 10) _phi45s = "0" + to_string(_phi45);	else _phi45s = to_string(_phi45);
		_phi46 = TMath::Max(int(floor(phi46/phistep)),0);
		if(_phi46 < 10) _phi46s = "0" + to_string(_phi46);	else _phi46s = to_string(_phi46);
		_phi56 = TMath::Max(int(floor(phi56/phistep)),0);
		if(_phi56 < 10) _phi56s = "0" + to_string(_phi56);	else _phi56s = to_string(_phi56);


		coords[i] = _m4s + _m5s + _m6s + _ue4s + _ue5s + _ue6s + _um4s + _um5s + _um6s + _phi45s + _phi46s + _phi56s + _chi2s;
		std::cout << float(i)/float(chi2_99->GetEntries()) << "\r";
	}
	std::cout << "\nCoords all full up.\nWe've got " <<  coords.size() << " of 'em" << std::endl;

	sort(coords.begin(),coords.end());
	coords.erase(unique(coords.begin(), coords.end(),stringCompare), coords.end());

	std::cout << "Holy shit, that worked.\nNow, fill the new ntuple" << std::endl;
	for(int i = 0; i < coords.size(); i++){
		// Gotta convert those boys from string coordinates into actual positions in phase space
		_m4 = stoi(coords[i].substr(0,2));
		float m4p = pow(10,(_m4/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		_m5 = stoi(coords[i].substr(2,2));
		float m5p = pow(10,(_m5/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
		_m6 = stoi(coords[i].substr(4,2));
		float m6p = pow(10,(_m6/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));

		_ue4 = stoi(coords[i].substr(6,2));
		float ue4p = _ue4/float(gridPoints)*(umax-umin);
		_ue5 = stoi(coords[i].substr(8,2));
		float ue5p = _ue5/float(gridPoints)*(umax-umin);
		_ue6 = stoi(coords[i].substr(10,2));
		float ue6p = _ue6/float(gridPoints)*(umax-umin);

		_um4 = stoi(coords[i].substr(12,2));
		float um4p = _um4/float(gridPoints)*(umax-umin);
		_um5 = stoi(coords[i].substr(14,2));
		float um5p = _um5/float(gridPoints)*(umax-umin);
		_um6 = stoi(coords[i].substr(16,2));
		float um6p = _um6/float(gridPoints)*(umax-umin);

		_phi45 = stoi(coords[i].substr(18,2));
		float phi45p = _phi45/float(gridPoints)*2*TMath::Pi();
		_phi46 = stoi(coords[i].substr(20,2));
		float phi46p = _phi46/float(gridPoints)*2*TMath::Pi();
		_phi56 = stoi(coords[i].substr(22,2));
		float phi56p = _phi56/float(gridPoints)*2*TMath::Pi();

		float chi2p = stof(coords[i].substr(24,6));

		if(steriles == 1)
			chi2_99_pr->Fill(chi2p,m4p,ue4p,um4p,0,0,0,0,0,0,0,0,0);
		if(steriles == 2)
			chi2_99_pr->Fill(chi2p,m4p,ue4p,um4p,m5p,ue5p,um5p,0,0,0,phi45p,0,0);
		if(steriles == 3)
			chi2_99_pr->Fill(chi2p,m4p,ue4p,um4p,m5p,ue5p,um5p,m6p,ue6p,um6p,phi45p,phi46p,phi56p);
	}
	std::cout << "99\% done!" << std::endl;

	// Now, for 90%, we don't need to do all that bullshit again. We can just go through 99% and take ones with appropriate chi2
	// find the chi2 min:
	chi2_99_pr->SetBranchAddress("chi2",&chi2);
	chi2_99_pr->SetBranchAddress("m4",&m4);
	chi2_99_pr->SetBranchAddress("ue4",&ue4);
	chi2_99_pr->SetBranchAddress("um4",&um4);
	chi2_99_pr->SetBranchAddress("m5",&m5);
	chi2_99_pr->SetBranchAddress("ue5",&ue5);
	chi2_99_pr->SetBranchAddress("um5",&um5);
	chi2_99_pr->SetBranchAddress("m6",&m6);
	chi2_99_pr->SetBranchAddress("ue6",&ue6);
	chi2_99_pr->SetBranchAddress("um6",&um6);
	chi2_99_pr->SetBranchAddress("phi45",&phi45);
	chi2_99_pr->SetBranchAddress("phi46",&phi46);
	chi2_99_pr->SetBranchAddress("phi56",&phi56);

	double chi2min = 3000;
	for(int i = 0; i < chi2_99->GetEntries(); i++){
	    chi2_99->GetEntry(i);
		if(chi2 < chi2min)
			chi2min = chi2;
	}
	std::cout << chi2min << std::endl;
	// Now fill up the 90%
	for(int i = 0; i < chi2_99_pr->GetEntries(); i++){
		chi2_99_pr->GetEntry(i);
		if(chi2 < chi2min + cl90)
			chi2_90_pr->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
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
	steriles = stoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line2(line);
	std::getline(is_line2,key,'=');
	std::getline(is_line2,value);
	nRuns = stoi(value.c_str());

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
