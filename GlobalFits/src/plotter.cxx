#include "TLegend.h"
#include "globalFit.h"
#include "TCut.h"
bool procOpt();

std::string plotOutput = "plots";
int steriles, nRuns, type, raster;
std::string dataset, location, output;
std::string procOptLoc;

int globFit_plotter(){

	procOptLoc = "/lar1nd/app/users/dcianci/SBN_3plusN/GlobalFits/inputs/";
    procOpt();

    std::cout << "Loading ntuple files..." << std::endl;
	std::string infile = output + Form("/nt_3%i_",steriles) + dataset + ".root";
	TString inputFile = infile;
	TFile *f = new TFile(inputFile);
	TNtuple *chi2_99 = (TNtuple*)(f->Get("chi2_99"));
	TNtuple *chi2_90 = (TNtuple*)(f->Get("chi2_90"));
	TNtuple *chi2_95 = (TNtuple*)(f->Get("chi2_95"));

    TCanvas *c1 = new TCanvas("c1");

	// Setup histo
    TH2F *h = new TH2F("h","3+3 #Chi^2;U_{e 4};U_{#mu 4}",1000,0.01,100.,1000,.01,100.);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(.8);
    h->GetXaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleFont(62);
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetLabelOffset(0.001);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.04);
    h->SetStats(kFALSE);

	c1->SetLogy();
	c1->SetLogx();
	chi2_99->SetMarkerStyle(7);
	chi2_90->SetMarkerStyle(7);
	chi2_99->SetFillColor(kBlue);
    chi2_90->SetFillColor(kMagenta);
	chi2_99->SetMarkerColor(kBlue);
	chi2_90->SetMarkerColor(kMagenta);
	chi2_95->SetMarkerStyle(7);
	chi2_95->SetMarkerColor(kRed+3);

	TLegend *leg = new TLegend(0.7,0.7,0.95,0.9);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextFont(62);
	leg->SetTextSize(0.03);
	leg->AddEntry(chi2_99,"99%% CL","f");
	leg->AddEntry(chi2_90,"90%% CL","f");

	if(steriles == 3){
		h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{51}");
		h->GetXaxis()->SetLimits(.01,100.);
		h->GetYaxis()->SetLimits(.01,100.);
		h->Draw();

		chi2_99->Draw("m5*m5:m4*m4","","same");
		chi2_90->Draw("m5*m5:m4*m4","","same");
		leg->Draw();
        c1->Print((plotOutput + "/" + dataset + "_dm251xdm241.png").c_str());

		h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{61}");
		h->GetXaxis()->SetLimits(.01,100.);
		h->GetYaxis()->SetLimits(.01,100.);
		h->Draw();

		chi2_99->Draw("(m6*m6):(m4*m4)","","same");
		chi2_90->Draw("(m6*m6):(m4*m4)","","same");
		leg->Draw();
        c1->Print(("plots/" + dataset + "_dm261xdm241.png").c_str());
	}

	if(steriles == 2){
		h->SetTitle("#chi^{2} for 3+2 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{51}");
		h->GetXaxis()->SetLimits(.01,100.);
		h->GetYaxis()->SetLimits(.01,100.);
		h->Draw();

		chi2_99->Draw("m5*m5:m4*m4","","same");
		chi2_90->Draw("m5*m5:m4*m4","","same");
		leg->Draw();
        c1->Print((plotOutput + "/" + dataset + "_3plus2_dm251xdm241.png").c_str());
	}

	if(steriles == 1){
		if(raster == 0){
			if(type==0)	h->SetTitle("#chi^{2} for 3+1 Sterile Fits;sin^{2}(2#Theta_{e#mu});#Delta m^{2}_{41}");
			if(type==1)	h->SetTitle("#chi^{2} for 3+1 Sterile Fits;sin^{2}(2#Theta_{#mu#mu});#Delta m^{2}_{41}");
			if(type==2)	h->SetTitle("#chi^{2} for 3+1 Sterile Fits;sin^{2}(2#Theta_{ee});#Delta m^{2}_{41}");
			h->GetXaxis()->SetLimits(.0001,.1);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			if(type==0){
        		chi2_99->Draw("m4*m4:4*um4*um4*um4*um4","","same");
        		chi2_90->Draw("m4*m4:4*um4*um4*um4*um4","","same");
			}
			else if(type == 1){
				chi2_99->Draw("m4*m4:4*um4*um4*(1-um4*um4)","","same");
        		chi2_90->Draw("m4*m4:4*um4*um4*(1-um4*um4)","","same");
			}
			else if(type == 2){
				chi2_99->Draw("m4*m4:4*ue4*ue4*(1-ue4*ue4)","","same");
        		chi2_90->Draw("m4*m4:4*ue4*ue4*(1-ue4*ue4)","","same");
			}
			leg->Draw();
        	c1->Print((plotOutput + "/" + dataset + "_3plus1_dm241xsinsq2t.png").c_str());
		}
		if(raster == 1){
			if(type==0)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{e#mu});#Delta m^{2}_{41}");
			if(type==1)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{#mu#mu});#Delta m^{2}_{41}");
			if(type==2)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{ee});#Delta m^{2}_{41}");
			h->GetXaxis()->SetLimits(.0001,.1);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_95->Draw("dm2:sin22th","","same");
			c1->Print((plotOutput + "/" + dataset + "_3plus1_dm241xsinsq2t.png").c_str());
		}
	}

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

    return true;
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
