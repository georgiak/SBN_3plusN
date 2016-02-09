#include "model.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"


SBN_spectrum::SBN_spectrum(struct neutrinoModel numodel){
	nullModel.zero();
	workingModel=numodel;	

	//Currently doesnt actually make a run. Just fudges for a difference
	//Will eventually auto oscillate!
	
	TFile *f = new TFile("data_hists/combined_ntuple_700m_numubar.root");
	TTree *t = (TTree*)f->Get("EventsTot");
	
	THmuboone_m = new TH1D("muboone_m","",N_m_bins,mu_bins);
	THsbnd_m = new TH1D("sbnd_m","",N_m_bins,mu_bins);
	THicarus_m = new TH1D("icarus_m","",N_m_bins,mu_bins);
	THmuboone_e = new TH1D("muboone_e","",N_e_bins,e_bins);
	THsbnd_e = new TH1D("sbnd_e","",N_e_bins,e_bins);
	THicarus_e = new TH1D("icarus_e","",N_e_bins,e_bins);	
	
	double energy;
	int parid;
	t->SetBranchAddress("energy",&energy);
	t->SetBranchAddress("parid",&parid);


	double FUDGE =0.9;

	for(int i=0; i<  t->GetEntries(); i++)
	{
		
		t->GetEntry(i);



		if(parid<20){
			THmuboone_m->Fill(energy,FUDGE);
			THsbnd_m->Fill(energy,FUDGE);
			THicarus_m->Fill(energy,FUDGE);
		}
		if(parid>20){
			THmuboone_e->Fill(energy,FUDGE);
			THsbnd_e->Fill(energy,FUDGE);
			THicarus_e->Fill(energy,FUDGE);
		}

	}

	//Dont actually oscillate yet, just smear!
/*	THmuboone_m
	THsbnd_m
	THicarus_m
	THmuboone_e
	THsbnd_e
	THicarus_e
*/



};

SBN_spectrum::SBN_spectrum(){
	nullModel.zero();
	workingModel=nullModel;
	
	TFile *f = new TFile("data_hists/combined_ntuple_700m_numubar.root");
	TTree *t = (TTree*)f->Get("EventsTot");
	
	THmuboone_m = new TH1D("muboone_m","",N_m_bins,mu_bins);
	THsbnd_m = new TH1D("sbnd_m","",N_m_bins,mu_bins);
	THicarus_m = new TH1D("icarus_m","",N_m_bins,mu_bins);
	THmuboone_e = new TH1D("muboone_e","",N_e_bins,e_bins);
	THsbnd_e = new TH1D("sbnd_e","",N_e_bins,e_bins);
	THicarus_e = new TH1D("icarus_e","",N_e_bins,e_bins);	
	
	double energy;
	int parid;
	t->SetBranchAddress("energy",&energy);
	t->SetBranchAddress("parid",&parid);


	for(int i=0; i<  t->GetEntries(); i++)
	{
		
		t->GetEntry(i);

		if(parid<20){
			THmuboone_m->Fill(energy);
			THsbnd_m->Fill(energy);
			THicarus_m->Fill(energy);
		}
		if(parid>20){
			THmuboone_e->Fill(energy);
			THsbnd_e->Fill(energy);
			THicarus_e->Fill(energy);
		}

	}
	//std::cout<<THmuboone_e->GetBinContent(0)<<std::endl;
};

void SBN_spectrum::THprint(){
	THmuboone_m->Print();
	THsbnd_m->Print();
	THicarus_m->Print();

	THmuboone_e->Print();
	THsbnd_e->Print();
	THicarus_e->Print();

//	TCanvas* c = new TCanvas();
	THmuboone_m->Print();
//	c->SaveAs("h.jpg");

};


const int SBN_spectrum::N_e_bins;
const int SBN_spectrum::N_m_bins;

constexpr double SBN_spectrum::mu_bins[N_m_bins+1];
constexpr double SBN_spectrum::e_bins[N_e_bins+1];

std::vector<double > SBN_spectrum::get_sixvector(){
			std::vector<double> ans = sbnd_e;
			ans.insert(std::end(ans), std::begin(muboone_e), std::end(muboone_e));
			ans.insert(std::end(ans), std::begin(icarus_e), std::end(icarus_e));
			ans.insert(std::end(ans), std::begin(sbnd_m), std::end(sbnd_m));
			ans.insert(std::end(ans), std::begin(muboone_m), std::end(muboone_m));
			ans.insert(std::end(ans), std::begin(icarus_m), std::end(icarus_m));
	
			return ans; 
};



void SBN_spectrum::update_model(struct neutrinoModel numodel){
	workingModel=numodel;
};

void SBN_spectrum::oscillate(){
//	muboone_e = BKG+ProbOscillate(MUBOONE_FLAG,E_FLAG);	
	//Will load the FULL_OSC, nu_E and nu_mu and eight them appropiately for each 
	//detector and each true L/E_v and so forth


};

void SBN_spectrum::fill_hists(){
	for(int i =0; i<muboone_e.size(); i++)
	{
		THmuboone_e->SetBinContent(i,muboone_e[i]);
		THsbnd_e->SetBinContent(i,sbnd_e[i]);
		THicarus_e->SetBinContent(i,icarus_e[i]);
	}

	for(int i =0; i<muboone_m.size(); i++)
	{
		THmuboone_m->SetBinContent(i,muboone_m[i]);
		THsbnd_m->SetBinContent(i,sbnd_m[i]);
		THicarus_m->SetBinContent(i,icarus_m[i]);
	}
};

void SBN_spectrum::fill_vectors(){
	for(int i =0; i<THmuboone_e->GetSize()-1; i++)
	{
		muboone_e.push_back( THmuboone_e->GetBinContent(i));
		sbnd_e.push_back(    THsbnd_e->GetBinContent(i));
		icarus_e.push_back(  THicarus_e->GetBinContent(i));
	}

	for(int i=0; i<THmuboone_m->GetSize()-1; i++)
	{
		muboone_m.push_back( THmuboone_m->GetBinContent(i));
		sbnd_m.push_back(    THsbnd_m->GetBinContent(i));
		icarus_m.push_back(  THicarus_m->GetBinContent(i));
	}
};

