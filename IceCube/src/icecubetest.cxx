#include "icecube.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int init(){

	//
	// First off: BRING IN THE FLUXES!
	//
	std::cout << "Initialize fluxes" << std::endl;

	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	double a_energy[700];
	double a_zenith[100];
	for(int i = 0; i < 100; i++)
		a_zenith[i] = -1. + i*1.2/100.;
	for(int i = 0; i < 700; i++){
		double logE = log10(1.) + i*(log10(1e6)-log10(1.))/700.;
		a_energy[i] = pow(10.0,logE);
	}
	double *a_flux = (double*)malloc(sizeof(a_energy) * sizeof(a_zenith) / sizeof(double));

	gsl_spline2d * spline_kaon = gsl_spline2d_alloc(T, 100, 700);
	gsl_spline2d * spline_pion = gsl_spline2d_alloc(T, 100, 700);
	gsl_interp_accel *accZen = gsl_interp_accel_alloc();
	gsl_interp_accel *accEn = gsl_interp_accel_alloc();

	// Load up our file
	ifstream file;
	file.open("libflux/flux_kaon_nominal.txt");
	// Skip first 12 lines because we don't need 'em!'
	file.ignore(10000,'\n');
	float dummy;
	double nufx, nubarfx;
	for(int cosz = 0; cosz < 100; cosz++){
    for(int recoe = 0; recoe < 700; recoe++){
			file >> dummy;
			file >> dummy;
			file >> nufx;
			file >> nubarfx;
			gsl_spline2d_set(spline_kaon, a_flux, cosz,recoe,nufx + nubarfx);
		}
	}
	file.close();

	file.open("libflux/flux_pion_nominal.txt");
	// Skip first 12 lines because we don't need 'em!'
	file.ignore(10000,'\n');
	for(int cosz = 0; cosz < 100; cosz++){
    for(int recoe = 0; recoe < 700; recoe++){
			file >> dummy;
			file >> dummy;
			file >> nufx;
			file >> nubarfx;
			gsl_spline2d_set(spline_pion, a_flux, cosz,recoe,nufx + nubarfx);
		}
	}
	file.close();

	gsl_spline2d_init(spline_kaon, a_zenith, a_energy, a_flux, 100, 700);
	gsl_spline2d_init(spline_pion, a_zenith, a_energy, a_flux, 100, 700);

	//
	// Next up: BRING IN AAAAALLLL THE MONTE CARLO AND WEIGHT IT TO OUR FLUX
	//
	std::cout << "Initializing MC" << std::endl;

	int nFOscEvts = 8784618;
	std::vector<float> mc_muonEnergyProxy, mc_recoMuonZenith, mc_weight;
	mc_muonEnergyProxy.resize(nFOscEvts);
	mc_recoMuonZenith.resize(nFOscEvts);
	mc_weight.resize(nFOscEvts);

	// Load up the nominal signal
	file.open("data/NuFSGenMC_nominal.dat");
	// Skip first 11 lines because we don't need 'em!'
	for(int i = 0; i < 11; i++)
		file.ignore(10000,'\n');

	// Read in the MC and weigh them according to our flux at icecube
	float weight, pion_flux, kaon_flux;
	for(int j = 0; j < nFOscEvts; j++){
		file >> dummy;
		file >> mc_muonEnergyProxy[j];
		file >> mc_recoMuonZenith[j];
		file >> dummy;
		file >> dummy;
		file >> weight;
		file >> dummy;
		file >> dummy;
		pion_flux = gsl_spline2d_eval(spline_pion,mc_recoMuonZenith[j],mc_muonEnergyProxy[j],accZen,accEn);
		kaon_flux = gsl_spline2d_eval(spline_kaon,mc_recoMuonZenith[j],mc_muonEnergyProxy[j],accZen,accEn);
		mc_weight[j] = (pion_flux+kaon_flux) * weight;
	}
	file.close();

	//
	// Now: TIME FOR DATAAAAAAAAA
	//
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

	//
	// Okay, we've made it this far. LET'S DRAW A HISTO AND SEE IF WE'RE CLOSE TO CORRECT
	//
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
	for(int i = 0; i < nEvts; i++){
		h->Fill(signal_muonEnergyProxy[i]);
	}
	for(int i = 0; i < nFOscEvts; i++){
		nom->Fill(mc_muonEnergyProxy[i],mc_weight[i]);
	}

	// Draw the histo to check that it all looks right.
	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetFillColor(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05);
	gStyle->SetHatchesLineWidth(2);
	gStyle->SetErrorX(0);

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
  h->SetMarkerStyle(4);

	c1->SetLogy();
	c1->SetLogx();

	h->Draw("PE");
	nom->Draw("H");
	c1->Print("IC8NomMC_Davio.pdf");

	gsl_spline2d_free(spline_kaon);
	gsl_spline2d_free(spline_pion);
	gsl_interp_accel_free(accZen);
	gsl_interp_accel_free(accEn);
	free(a_flux);

	return 0;
}



int main()
{
	init();
    return 0;
}
