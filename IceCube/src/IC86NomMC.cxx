#include "icecube.h"

#include <gsl/gsl_math.h>

int init(){

	// Load up our file
	int nFOscEvts = 8784618;
	std::vector<float> mc_muonEnergyProxy, mc_recoMuonZenith, mc_weight;
	mc_muonEnergyProxy.resize(nFOscEvts);
	mc_recoMuonZenith.resize(nFOscEvts);
	mc_weight.resize(nFOscEvts);

	// Load up the nominal MC
  std::cout << "Initialize MC" << std::endl;
  ifstream file;
	file.open("data/NuFSGenMC_nominal.dat");
	// Skip first 11 lines because we don't need 'em!'
	for(int i = 0; i < 11; i++)
		file.ignore(10000,'\n');

	// Read in the MC and weigh them according to our flux at icecube
	float weight, pion_flux, kaon_flux, dummy;
	for(int j = 0; j < nFOscEvts; j++){
		file >> dummy;
		file >> mc_muonEnergyProxy[j];
		file >> mc_recoMuonZenith[j];
		file >> dummy;
		file >> dummy;
		file >> weight;
		file >> pion_flux;
		file >> kaon_flux;
		mc_weight[j] = (pion_flux+kaon_flux) * weight;
	}
	file.close();


  // Load up the data
	std::cout << "Initialize Data" << std::endl;

	int nEvts = 20145;
	std::vector<float> signal_muonEnergyProxy, signal_recoMuonZenith;
	signal_muonEnergyProxy.resize(nEvts);
	signal_recoMuonZenith.resize(nEvts);

	// Load up the data
	file.open("data/observed_events.dat");
	// Skip first 12 lines because we don't need 'em!'
	for(int i = 0; i < 12; i++)
		file.ignore(10000,'\n');

	for(int j = 0; j < nEvts; j++){
		file >> signal_muonEnergyProxy[j];
		file >> signal_recoMuonZenith[j];
	}
	file.close();


  // Draw plots to see how we're doing
	std::cout << "Get plottin'" << std::endl;

	//I want evenly spaced log bins
	float emin = 400;	float emax = 20000;
	int nbins = 28;
	float bins[nbins+1];
	for(int i = 0; i < nbins+1; i++){
		bins[i] = pow(10,log10(emin) + i*log10(emax/emin)/nbins);
	}

	TH1F *h = new TH1F("h","Reconstructed Events;Energy/GeV",nbins,bins);
	TH1F *nom = new TH1F("nom","nom",nbins,bins);
  TH1F *ratio = new TH1F("ratio","ratio",nbins,bins);
	for(int i = 0; i < nEvts; i++){
		h->Fill(signal_muonEnergyProxy[i]);
	}
	for(int i = 0; i < nFOscEvts; i++){
		nom->Fill(mc_muonEnergyProxy[i],mc_weight[i]);
	}
  float sig,mc;
  for(int b = 1; b < nbins+1; b++){
    sig = h->GetBinContent(b);
    mc = nom->GetBinContent(b);
    ratio->SetBinContent(b,sig/mc);
    //std::cout << "BIN: " << b << " RATIO: " << sig/mc << std::endl;
  }

/*
	// Now, we've gotta normalize the bins of the nominal distribution to match the signal.
	int datact = 0;	int nomct = 0;
	for(int b = 1; b < nbins+1; b++){
		datact += h->GetBinContent(b);
		nomct += nom->GetBinContent(b);
	}
	for(int i = 1; i < nbins+1; i++){
		nom->SetBinContent(i,nom->GetBinContent(i)*datact/nomct);
	}
*/

	// Draw the histo to check that it all looks right.
	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetFillColor(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05);
  gStyle->SetErrorX(0);
	gStyle->SetHatchesLineWidth(2);

	h->GetXaxis()->SetTitleOffset(1.1);
	h->GetYaxis()->SetTitleOffset(.8);
	h->GetXaxis()->SetTitleFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->CenterTitle();
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleSize(0.04);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetXaxis()->SetLabelOffset(0.001);
	h->GetYaxis()->SetTitleSize(0.04);
	h->GetYaxis()->SetLabelSize(0.04);
	h->SetStats(kFALSE);
	h->SetMarkerColor(2);
  h->SetMarkerStyle(kFullCircle);
  nom->SetLineColor(kBlack);

	c1->SetLogy();
	c1->SetLogx();

	h->Draw("PE");
	nom->Draw("H same");
	c1->Print("IC86NomMC.pdf");



  // Draw Ratio Plot for test
  TCanvas *c2 = new TCanvas("c2");
  ratio->GetXaxis()->SetTitleOffset(1.1);
	ratio->GetYaxis()->SetTitleOffset(.8);
	ratio->GetXaxis()->SetTitleFont(62);
	ratio->GetYaxis()->SetTitleFont(62);
	ratio->GetYaxis()->CenterTitle();
	ratio->GetXaxis()->CenterTitle();
	ratio->GetXaxis()->SetTitleSize(0.04);
	ratio->GetXaxis()->SetLabelSize(0.04);
	ratio->GetXaxis()->SetLabelOffset(0.001);
	ratio->GetYaxis()->SetTitleSize(0.04);
	ratio->GetYaxis()->SetLabelSize(0.04);
	ratio->SetStats(kFALSE);
	ratio->SetMarkerColor(2);
  ratio->SetMarkerStyle(kFullCircle);
  ratio->Draw("P");
  c2->SetLogx();
  c2->Print("IC98NomMC_ratio.pdf");


	return 0;
}

int main()
{
	init();
    return 0;
}
