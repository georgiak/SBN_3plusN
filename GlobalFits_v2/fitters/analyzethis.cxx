#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH3F.h"
#include "globalFit.h"
#include "TCut.h"
#include "TView3D.h"
#include "TGraphPainter.h"
bool procOpt();

std::string plotOutput = "plots";
int steriles, nRuns, type, raster, discretized, diag, dims;
std::string dataset, location, output;
std::string procOptLoc;
std::string suffix;
float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,step,temp;
float m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;

float onesigma[] = {3.53, 8.18, 13.74};

int globFit_analyze(){

	procOptLoc = "inputs/";
    procOpt();

	TMarker *bestfit = new TMarker();
	TH1D *h_chi2 = new TH1D("chi2","chi2;chi2",1000,200,300);
	TH1D *h_step = new TH1D("step","step;step",100,0,.2);
	TH1D *h_temp = new TH1D("temp","temp;temp",100,0,2);
	TH1D *h_m4 = new TH1D("m4","m4;eV",100,.1,10);
	TH1D *h_ue4 = new TH1D("ue4","ue4;ue",100,.01,.5);
	TH1D *h_um4 = new TH1D("um4","um4;um",100,.01,.5);
	TH1D *h_m5 = new TH1D("m5","m5;eV",100,.1,10);
	TH1D *h_ue5 = new TH1D("ue5","ue5;ue",100,.01,.5);
	TH1D *h_um5 = new TH1D("um5","um5;um",100,.01,.5);
	TH1D *h_m6 = new TH1D("m6","m6;eV",100,.1,10);
	TH1D *h_ue6 = new TH1D("ue6","ue6;ue",100,.01,.5);
	TH1D *h_um6 = new TH1D("um6","um6;um",100,.01,.5);
	TH1D *h_phi45 = new TH1D("phi45","phi45;Radians",100,0,2*TMath::Pi());
	TH1D *h_phi46 = new TH1D("phi46","phi46;Radians",100,0,2*TMath::Pi());
	TH1D *h_phi56 = new TH1D("phi56","phi56;Radians",100,0,2*TMath::Pi());

    std::cout << "Loading ntuple files..." << std::endl;
	std::string infile;
	if(discretized == 0)	infile = output + Form("/nt_3%i_",steriles) + dataset + ".root";
	if(discretized == 1)	infile = output + Form("/nt_3%i_",steriles) + dataset + "_processed.root";
	TString inputFile = infile;
	std::cout << "Infile: " << infile << std::endl;
	TFile *f = new TFile(inputFile);
	TNtuple *chi2_99_all;
	TNtuple *chi2_90_all;
	TNtuple *chi2_95;

	float um4min, um4max, um5min, um5max, um6min, um6max, ue4min, ue4max, ue5min, ue5max, ue6min, ue6max;
	float dm41min, dm41max, dm51min, dm51max, dm61min, dm61max;
	um4min = 1.; um4max = 0.;
	um5min = 1.; um5max = 0.;
	um6min = 1.; um6max = 0.;
	ue4min = 1.; ue4max = 0.;
	ue5min = 1.; ue5max = 0.;
	ue6min = 1.; ue6max = 0.;
	dm41min = 1.;	dm41max = 0.;
	dm51min = 1.;	dm51max = 0.;
	dm61min = 1.;	dm61max = 0.;

	if(discretized == 0){
		chi2_99_all = (TNtuple*)(f->Get("chi2_99"));
		chi2_90_all = (TNtuple*)(f->Get("chi2_90"));
		chi2_95 = (TNtuple*)(f->Get("chi2_95"));
		suffix = "";
	}
	if(discretized == 1){
		chi2_99_all = (TNtuple*)(f->Get("chi2_99_pr"));
		chi2_90_all = (TNtuple*)(f->Get("chi2_90_pr"));
		suffix = "_disc";
	}

	// Find chi2Min
	chi2_99_all->SetBranchAddress("chi2",&chi2);
	if(dims == 1){
		chi2_99_all->SetBranchAddress("step",&step);
		chi2_99_all->SetBranchAddress("temp",&temp);
	}
	chi2_99_all->SetBranchAddress("m4",&m4);
	chi2_99_all->SetBranchAddress("ue4",&ue4);
	chi2_99_all->SetBranchAddress("um4",&um4);
	chi2_99_all->SetBranchAddress("m5",&m5);
	chi2_99_all->SetBranchAddress("ue5",&ue5);
	chi2_99_all->SetBranchAddress("um5",&um5);
	chi2_99_all->SetBranchAddress("m6",&m6);
	chi2_99_all->SetBranchAddress("ue6",&ue6);
	chi2_99_all->SetBranchAddress("um6",&um6);
	chi2_99_all->SetBranchAddress("phi45",&phi45);
	chi2_99_all->SetBranchAddress("phi46",&phi46);
	chi2_99_all->SetBranchAddress("phi56",&phi56);
	float chi2min = 3000.f;
    for(int i = 0; i < chi2_99_all->GetEntries(); i++){
        chi2_99_all->GetEntry(i);
		if(chi2 < chi2min){
			chi2min = chi2;
        	m4_min = m4;	ue4_min = ue4;	um4_min = um4;
			m5_min = m5;	ue5_min = ue5;	um5_min = um5;
			m6_min = m6;	ue6_min = ue6;	um6_min = um6;
			phi45_min = phi45;	phi46_min = phi46;	phi56_min = phi56;
		}
		um4min = min(um4min, um4); um4max = max(um4max, um4);
		um5min = min(um5min, um5); um5max = max(um5max, um5);
		um6min = min(um6min, um6); um6max = max(um6max, um6);
		ue4min = min(ue4min, ue4); ue4max = max(ue4max, ue4);
		ue5min = min(ue5min, ue5); ue5max = max(ue5max, ue5);
		ue6min = min(ue6min, ue6); ue6max = max(ue6max, ue6);
		dm41min = min(dm41min, m4); dm41max = max(dm41max, m4);
		dm51min = min(dm51min, m5); dm51max = max(dm51max, m5);
		dm61min = min(dm61min, m6); dm61max = max(dm61max, m6);
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
	std::cout << "For 99\%" << std::endl;
	std::cout << "um4min: " << um4min << " -- um4max: " << um4max << std::endl;
	std::cout << "ue4min: " << ue4min << " -- ue4max: " << ue4max << std::endl;
	std::cout << "um5min: " << um5min << " -- um5max: " << um5max << std::endl;
	std::cout << "ue5min: " << ue5min << " -- ue5max: " << ue5max << std::endl;
	std::cout << "um6min: " << um6min << " -- um6max: " << um6max << std::endl;
	std::cout << "ue6min: " << ue6min << " -- ue6max: " << ue6max << std::endl;
	std::cout << "dm241min: " << pow(dm41min,2) << " -- dm241max: " << pow(dm41max,2) << std::endl;
	std::cout << "dm251min: " << pow(dm51min,2) << " -- dm251max: " << pow(dm51max,2) << std::endl;
	std::cout << "dm261min: " << pow(dm61min,2) << " -- dm261max: " << pow(dm61max,2) << std::endl;

	chi2_90_all->SetBranchAddress("chi2",&chi2);
	chi2_90_all->SetBranchAddress("m4",&m4);
	chi2_90_all->SetBranchAddress("ue4",&ue4);
	chi2_90_all->SetBranchAddress("um4",&um4);
	chi2_90_all->SetBranchAddress("m5",&m5);
	chi2_90_all->SetBranchAddress("ue5",&ue5);
	chi2_90_all->SetBranchAddress("um5",&um5);
	chi2_90_all->SetBranchAddress("m6",&m6);
	chi2_90_all->SetBranchAddress("ue6",&ue6);
	chi2_90_all->SetBranchAddress("um6",&um6);
	chi2_90_all->SetBranchAddress("phi45",&phi45);
	chi2_90_all->SetBranchAddress("phi46",&phi46);
	chi2_90_all->SetBranchAddress("phi56",&phi56);
	um4min = 1.; um4max = 0.;
	um5min = 1.; um5max = 0.;
	um6min = 1.; um6max = 0.;
	ue4min = 1.; ue4max = 0.;
	ue5min = 1.; ue5max = 0.;
	ue6min = 1.; ue6max = 0.;
	dm41min = 1.;	dm41max = 0.;
	dm51min = 1.;	dm51max = 0.;
	dm61min = 1.;	dm61max = 0.;
	for(int i = 0; i < chi2_90_all->GetEntries(); i++){
        chi2_90_all->GetEntry(i);
		um4min = min(um4min, um4); um4max = max(um4max, um4);
		um5min = min(um5min, um5); um5max = max(um5max, um5);
		um6min = min(um6min, um6); um6max = max(um6max, um6);
		ue4min = min(ue4min, ue4); ue4max = max(ue4max, ue4);
		ue5min = min(ue5min, ue5); ue5max = max(ue5max, ue5);
		ue6min = min(ue6min, ue6); ue6max = max(ue6max, ue6);
		dm41min = min(dm41min, m4); dm41max = max(dm41max, m4);
		dm51min = min(dm51min, m5); dm51max = max(dm51max, m5);
		dm61min = min(dm61min, m6); dm61max = max(dm61max, m6);
    }
	std::cout << "For 90\%" << std::endl;
	std::cout << "um4min: " << um4min << " -- um4max: " << um4max << std::endl;
	std::cout << "ue4min: " << ue4min << " -- ue4max: " << ue4max << std::endl;
	std::cout << "um5min: " << um5min << " -- um5max: " << um5max << std::endl;
	std::cout << "ue5min: " << ue5min << " -- ue5max: " << ue5max << std::endl;
	std::cout << "um6min: " << um6min << " -- um6max: " << um6max << std::endl;
	std::cout << "ue6min: " << ue6min << " -- ue6max: " << ue6max << std::endl;
	std::cout << "dm241min: " << pow(dm41min,2) << " -- dm241max: " << pow(dm41max,2) << std::endl;
	std::cout << "dm251min: " << pow(dm51min,2) << " -- dm251max: " << pow(dm51max,2) << std::endl;
	std::cout << "dm261min: " << pow(dm61min,2) << " -- dm261max: " << pow(dm61max,2) << std::endl;

	// Now do 1sigma
	um4min = 1.; um4max = 0.;
	um5min = 1.; um5max = 0.;
	um6min = 1.; um6max = 0.;
	ue4min = 1.; ue4max = 0.;
	ue5min = 1.; ue5max = 0.;
	ue6min = 1.; ue6max = 0.;
	dm41min = 1.;	dm41max = 0.;
	dm51min = 1.;	dm51max = 0.;
	dm61min = 1.;	dm61max = 0.;
	for(int i = 0; i < chi2_90_all->GetEntries(); i++){
        chi2_90_all->GetEntry(i);
		if(chi2 - chi2min < onesigma[steriles-1]){
			um4min = min(um4min, um4); um4max = max(um4max, um4);
			um5min = min(um5min, um5); um5max = max(um5max, um5);
			um6min = min(um6min, um6); um6max = max(um6max, um6);
			ue4min = min(ue4min, ue4); ue4max = max(ue4max, ue4);
			ue5min = min(ue5min, ue5); ue5max = max(ue5max, ue5);
			ue6min = min(ue6min, ue6); ue6max = max(ue6max, ue6);
			dm41min = min(dm41min, m4); dm41max = max(dm41max, m4);
			dm51min = min(dm51min, m5); dm51max = max(dm51max, m5);
			dm61min = min(dm61min, m6); dm61max = max(dm61max, m6);
		}
    }
	std::cout << "For 1sigma" << std::endl;
	std::cout << "um4min: " << um4min << " -- um4max: " << um4max << std::endl;
	std::cout << "ue4min: " << ue4min << " -- ue4max: " << ue4max << std::endl;
	std::cout << "um5min: " << um5min << " -- um5max: " << um5max << std::endl;
	std::cout << "ue5min: " << ue5min << " -- ue5max: " << ue5max << std::endl;
	std::cout << "um6min: " << um6min << " -- um6max: " << um6max << std::endl;
	std::cout << "ue6min: " << ue6min << " -- ue6max: " << ue6max << std::endl;
	std::cout << "dm241min: " << pow(dm41min,2) << " -- dm241max: " << pow(dm41max,2) << std::endl;
	std::cout << "dm251min: " << pow(dm51min,2) << " -- dm251max: " << pow(dm51max,2) << std::endl;
	std::cout << "dm261min: " << pow(dm61min,2) << " -- dm261max: " << pow(dm61max,2) << std::endl;


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

	std::getline(file,line);
	std::istringstream is_line6(line);
	std::getline(is_line6,key,'=');
	std::getline(is_line6,value);
	type = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line7(line);
	std::getline(is_line7,key,'=');
	std::getline(is_line7,value);
	raster = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line8(line);
	std::getline(is_line8,key,'=');
	std::getline(is_line8,value);
	discretized = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line9(line);
	std::getline(is_line9,key,'=');
	std::getline(is_line9,value);
 	diag = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line10(line);
	std::getline(is_line10,key,'=');
	std::getline(is_line10,value);
 	dims = atoi(value.c_str());

    return true;
}

#ifndef __CINT__
int main()
{
    globFit_analyze();
    return 0;
}
# endif

void plotter(){
    globFit_analyze();
    return;
}
