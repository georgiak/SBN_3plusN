#include "TLegend.h"
#include "globalFit.h"
#include "TCut.h"

int steriles = 3;
std::string dataset = "all";
std::string myRoot = "nt_31_all.root";
std::string location = "ntuples";
std::string output = "plots";

int globFit_plotter(){

    std::cout << "Loading ntuple files..." << std::endl;
	std::string infile = location + "/" + myRoot;
	TString inputFile = infile;
	TFile *f = new TFile(inputFile);
	TNtuple *chi2_99 = (TNtuple*)(f->Get("chi2_99"));
	TNtuple *chi2_90 = (TNtuple*)(f->Get("chi2_90"));

    TCanvas *c1 = new TCanvas("c1");

	// Setup histo
    TH2F *h = new TH2F("h","3+3 #Chi^2;U_{e 4};U_{#mu 4}",1000,0.01,100.,1000,.01,100.);
    h->GetXaxis()->SetTitleOffset(.8);
    h->GetYaxis()->SetTitleOffset(.8);
    h->GetXaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetLabelOffset(0.001);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.04);
    h->SetStats(kFALSE);

	if(steriles == 3){
		c1->SetLogy();
		c1->SetLogx();

		h->SetTitle("3+3 #Chi^{2};#Delta m^{2}_{41};#Delta m^{2}_{51}");
		//h->GetXaxis()->SetRange(.01,100.);
		//h->GetYaxis()->SetRange(.01,100.);
		h->Draw();
        c1->Update();
        chi2_99->SetMarkerStyle(20);
		chi2_90->SetMarkerStyle(20);
        chi2_99->SetMarkerColor(kBlue);    		chi2_99->Draw("m5*m5:m4*m4","","same");
		chi2_90->SetMarkerColor(kMagenta);    	chi2_90->Draw("m5*m5:m4*m4","","same");
        c1->Print((output + "/" + dataset + "_dm251xdm241.png").c_str());

		h->SetTitle("3+3 #Chi^{2};#Delta m^{2}_{41};#Delta m^{2}_{61}");
		//h->GetXaxis()->SetRange(.01,100.);
		//h->GetYaxis()->SetRange(.01,100.);
		h->Draw();
        c1->Update();
		chi2_99->SetMarkerStyle(20);
		chi2_90->SetMarkerStyle(20);
        chi2_99->SetMarkerColor(kBlue);    		chi2_99->Draw("(m6*m6):(m4*m4)","","same");
		chi2_90->SetMarkerColor(kMagenta);    	chi2_90->Draw("(m6*m6):(m4*m4)","","same");
        c1->Print(("plots/" + dataset + "_dm261xdm241.png").c_str());
	}

	if(steriles == 2){
		c1->SetLogy();
		c1->SetLogx();

		h->SetTitle("3+2 #Chi^{2};#Delta m^{2}_{41};#Delta m^{2}_{51}");
		//h->GetXaxis()->SetRange(.01,100.);
		//h->GetYaxis()->SetRange(.01,100.);
		h->Draw();
        c1->Update();
		chi2_99->SetMarkerStyle(20);
		chi2_90->SetMarkerStyle(20);
        chi2_99->SetMarkerColor(kBlue);    		chi2_99->Draw("m5*m5:m4*m4","","same");
		chi2_90->SetMarkerColor(kMagenta);    	chi2_90->Draw("m5*m5:m4*m4","","same");
        c1->Print((output + "/" + dataset + "_3plus2_dm251xdm241.png").c_str());
	}

	if(steriles == 1){
		c1->SetLogy();

		h->SetTitle("3+1 #Chi^{2};sin^{2}(2#Theta_{e#mu});#Delta m^{2}_{41}")
		h->GetXaxis()->SetRange(.01,100.);
		h->GetYaxis()->SetRange(.01,100.);
		h->Draw();
        c1->Update();
        in_chain->SetMarkerStyle(20);
        in_chain->SetMarkerColor(kBlue);    	in_chain->Draw("m4*m4:4*ue4*ue4*um4*um4",cutS1.c_str(),"same");
        in_chain->SetMarkerColor(kMagenta);    	in_chain->Draw("m4*m4:4*ue4*ue4*um4*um4",cutS2.c_str(),"same");
        c1->Print(output + "/" + dataset + "_3plus1_dm241xsinsq2t.png");
	}


    return 0;
}


#ifndef __CINT__
int main()
{
    globFit_plotter();
    return 0;
}
# endif

void plotter(){
    globFit_plotter();
    return;
}
