#include "TF1.h"
#include "/home/mark/programs/root/root/include/Math/GSLIntegrator.h"
#include "/home/mark/programs/root/root/include/Math/WrappedTF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "/home/mark/programs/root/root/include/TRandom.h"
#include <iostream>
#include "THStack.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TLegend.h"
#include <vector>
#include "/home/mark/programs/root/root/include/Math/GSLIntegrator.h"
#include "/home/mark/programs/root/root/include/Math/WrappedTF1.h"


void plot_uboone_intrinsic(){


	double UBOONE_proposal_intrinsic[11] ={237.789, 428.737, 523.626, 554.084, 558.77, 509.569, 450.995, 370.164, 297.533, 215.531, 107.756};
	double UBOONE_proposal_photon[11] ={194.463, 125.347, 59.7447, 32.801, 19.9149, 14.0576, 10.5432, 9.37172, 8.20025, 9.37172, 2.34293};
	double UBOONE_proposal_muon[11] ={42.1727, 43.3442, 25.7722, 22.2578, 18.7434, 15.229, 11.7146, 10.5432, 7.02879, 4.68586, 4.68586};

	static const int N_m_bins = 20;
	static const int N_e_bins = 11;
	static  double mu_bins[N_m_bins+1] = {0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.25,1.5,2,2.5,3};
	static  double e_bins[N_e_bins+1]= {0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.3,1.5,1.75,2,3};	 
	static  double ebinw[N_e_bins]   = {0.15,0.15,0.15,0.15,0.15,0.15,0.2,0.2,0.25,0.25,1.0};
	static  double mubinw[N_m_bins]  = {0.1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.25,0.25,0.5,0.5,0.5};

	TH1D * Residual_intrinsic = new TH1D("Residual_intrinsic","",11,e_bins);
	TH1D * Residual_photon = new TH1D("Residual_photon","",11,e_bins);
	TH1D * Residual_muon = new TH1D("Residual_muon","",11,e_bins);
	
	TFile *f = new TFile("../bkg_data/uBooNE_intrin.root"); 
 
	TH1D * H1nue_reco_intrinsic = (TH1D*)f->Get("detector_intrinsic_app");
	TH1D * H1nue_reco_muon= (TH1D*)f->Get("detector_muon_app");
	TH1D * H1nue_reco_photon = (TH1D*)f->Get("detector_photon_app");

	gStyle->SetOptStat(0);

	double Nintrin = H1nue_reco_intrinsic->GetSumOfWeights();
	double Nmuon = H1nue_reco_muon->GetSumOfWeights();
	double Nphoton= H1nue_reco_photon->GetSumOfWeights();
	double ti =0;
	double ti2=0;
	double tp =0;
	double tp2 = 0;
	double tm = 0;
	double tm2 = 0;
	int k =0;


	for(k=0; k<11;k++)
	{
		ti = H1nue_reco_intrinsic->GetBinContent(k+1);
		tm = H1nue_reco_muon->GetBinContent(k+1);
		tp = H1nue_reco_photon->GetBinContent(k+1);

		Nintrin += ti;
		Nmuon += tm;
		Nphoton += tp;

	    /*     ti2 =100*(ti-UBOONE_proposal_intrinsic[k]*ebinw[k])/(338.0+396+94);
		 tp2 =100*(tp-UBOONE_proposal_photon[k]*ebinw[k])/(81+5.0);
		 tm2 =100*(tm-UBOONE_proposal_muon[k]*ebinw[k])/(35);
	*/
	     	ti2 =100*(ti-UBOONE_proposal_intrinsic[k]*ebinw[k])/(UBOONE_proposal_intrinsic[k]*ebinw[k]);
		 tp2 =100*(tp-UBOONE_proposal_photon[k]*ebinw[k])/(UBOONE_proposal_photon[k]*ebinw[k]);
		 tm2 =100*(tm-UBOONE_proposal_muon[k]*ebinw[k])/(UBOONE_proposal_muon[k]*ebinw[k]);

		Residual_intrinsic->SetBinContent(k+1,ti2);
		Residual_photon->SetBinContent(k+1,tp2);
		Residual_muon->SetBinContent(k+1,tm2);

		H1nue_reco_intrinsic->SetBinContent(k+1,ti/ebinw[k]);
		H1nue_reco_photon->SetBinContent(k+1,tp/ebinw[k]);
		H1nue_reco_muon->SetBinContent(k+1,tm/ebinw[k]);

	}

char namei[200];
sprintf(namei, "CC intrinsic #nu_{e} percen: %.1f ",100*Nintrin/(338.0+396+94));
	char namep[200];
sprintf(namep, "NC single #gamma percen: %.1f ",100*Nphoton/86.0);
	char namem[200];
sprintf(namem, "CC #mu mis-id percen: %.1f ",100*Nmuon/35.0);





	std::cout<<"UBOONE analysis"<<std::endl;
	std::cout<<"Nintrinsic: "<<Nintrin<<" proposal: 800ish"<<std::endl;
	std::cout<<"Nmuon: "<<Nmuon<<" proposal: 35"<<std::endl;
	std::cout<<"Nphoton: "<<Nphoton<<" proposal: 86"<<std::endl;


	THStack *hs = new THStack("hs","Stacked 1D histograms");
	H1nue_reco_intrinsic->SetFillColor(kGreen-6);
	H1nue_reco_intrinsic->SetMarkerStyle(21);
	H1nue_reco_intrinsic->SetMarkerColor(kGreen-6);
	H1nue_reco_intrinsic->SetLineColor(kBlack);


	H1nue_reco_photon->SetFillColor(kRed-7);
   	H1nue_reco_photon->SetMarkerStyle(21);
      	H1nue_reco_photon->SetMarkerColor(kRed-7);
      	H1nue_reco_photon->SetLineColor(kBlack);
	H1nue_reco_muon->SetFillColor(kBlue-7);
   	H1nue_reco_muon->SetMarkerStyle(21);
      	H1nue_reco_muon->SetMarkerColor(kBlue-7);	
      	H1nue_reco_muon->SetLineColor(kBlack);

	hs->Add(H1nue_reco_intrinsic);
	hs->Add(H1nue_reco_photon);
	hs->Add(H1nue_reco_muon);


	TCanvas* Cstack = new TCanvas();
	Cstack->SetCanvasSize(500*1.2,900*1.2);
	Cstack->Divide(1,3);
	Cstack->cd(1);
		hs->SetTitle("UBOONE, 13.2e20 POT ");
		hs->Draw();

		hs->GetXaxis()->SetTitle("Recontructed Energy (GeV)");
		hs->GetYaxis()->SetTitle("Events / GeV");

		hs->GetYaxis()->SetTitleOffset(1.5);

//		hs->GetYaxis()->SetRangeUser(0, 20000);
		hs->SetMaximum(1000);
		Cstack->Update();

		TLegend * legStack = new TLegend(0.6,0.65,0.85,0.94);
		legStack->AddEntry(H1nue_reco_intrinsic,namei,"f");
		legStack->AddEntry(H1nue_reco_photon,namep,"f");
		legStack->AddEntry(H1nue_reco_muon,namem,"f");
		legStack->Draw();	


	Cstack->cd(3);
	
	TLine *l=new TLine(0.2,0,3.0,0);
	l->SetLineColor(kBlack);
Residual_intrinsic->SetMaximum(15);
Residual_intrinsic->SetMinimum(-15);
	Residual_intrinsic->SetLineColor(kGreen-6);
	Residual_photon->SetLineColor(kRed-7);
	Residual_muon->SetLineColor(kBlue-7);

	Residual_muon->SetLineWidth(2);
	Residual_photon->SetLineWidth(2);
	Residual_intrinsic->SetLineWidth(2);
	Residual_intrinsic->SetTitle("Percentage Residual of SBN Proposal");
	Residual_intrinsic->GetXaxis()->SetTitle("Recontructed Energy (GeV)");
	Residual_intrinsic->GetYaxis()->SetTitle("% Residual");

	Residual_intrinsic->Draw();
	l->Draw("same");

	Cstack->cd(2);
	Residual_photon->SetTitle("Percentage Residual of SBN Proposal");
	Residual_photon->GetXaxis()->SetTitle("Recontructed Energy (GeV)");
	Residual_photon->GetYaxis()->SetTitle("% Residual");
Residual_photon->SetMaximum(200);
Residual_photon->SetMinimum(-200);

	Residual_photon->Draw();
	Residual_muon->Draw("same");

	l->Draw("same");

		Cstack->SaveAs("uBooNE_reco_stack.png");
	
	
}



/***********************TO DO***
 *
 *
 *make plots pretty and off one color
 *muon and pion decays, should be easy.
 *maybe a nice thing showing positional hits?
 *
 *MUboone, uBooNE and icarus all taken care of.
 *
 *
 *
 *
 *
 *
 *
 * */
