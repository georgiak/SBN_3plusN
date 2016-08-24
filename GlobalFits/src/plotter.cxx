#include "TLegend.h"
#include "TH3F.h"
#include "globalFit.h"
#include "TCut.h"
#include "TView3D.h"
bool procOpt();

std::string plotOutput = "plots";
int steriles, nRuns, type, raster, discretized, diag, dims;
std::string dataset, location, output;
std::string procOptLoc;
std::string suffix;
float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
float m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;

bool twodof = false;

int globFit_plotter(){

	procOptLoc = "/Users/dcianci/Physics/SBN_3plusN/GlobalFits/inputs/";
    procOpt();

	TH1D *h_chi2 = new TH1D("chi2","chi2;chi2",1000,200,300);
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

	if(discretized == 0){
		chi2_99_all = (TNtuple*)(f->Get("chi2_99"));
		chi2_90_all = (TNtuple*)(f->Get("chi2_90"));
		suffix = "";
	}
	if(discretized == 1){
		chi2_99_all = (TNtuple*)(f->Get("chi2_99_pr"));
		chi2_90_all = (TNtuple*)(f->Get("chi2_90_pr"));
		suffix = "_disc";
	}

	// Find chi2Min
	chi2_99_all->SetBranchAddress("chi2",&chi2);
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
		h_chi2->Fill(chi2);
		h_m4->Fill(m4);			h_ue4->Fill(ue4);		h_um4->Fill(um4);
		h_m5->Fill(m5);			h_ue5->Fill(ue5);		h_um5->Fill(um5);
		h_m6->Fill(m6);			h_ue6->Fill(ue6);		h_um6->Fill(um6);
		h_phi45->Fill(phi45);		h_phi46->Fill(phi46);		h_phi56->Fill(phi56);
		if(chi2 < chi2min){
			chi2min = chi2;
        	m4_min = m4;	ue4_min = ue4;	um4_min = um4;
			m5_min = m5;	ue5_min = ue5;	um5_min = um5;
			m6_min = m6;	ue6_min = ue6;	um6_min = um6;
			phi45_min = phi45;	phi46_min = phi46;	phi56_min = phi56;
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

	// Now, if we want, we can switch to only two dof for the chi2, which is what we want for these plots
	TNtuple *chi2_99 = new TNtuple("chi2Nt","chi2Nt","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	TNtuple *chi2_90 = new TNtuple("chi2Nt","chi2Nt","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	for(int i = 0; i < chi2_99_all->GetEntries(); i++){
        chi2_99_all->GetEntry(i);
		if(chi2-chi2min < 9.21){
			chi2_99->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
		}
		if(chi2-chi2min < 4.61){
			chi2_90->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
		}
	}
	std::cout << chi2_99->GetEntries() << " " << chi2_99_all->GetEntries() << std::endl;
	std::cout << chi2_90->GetEntries() << " " << chi2_90_all->GetEntries() << std::endl;


	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetFillColor(0);
  	gStyle->SetPadLeftMargin(0.15);
  	gStyle->SetPadBottomMargin(0.15);
  	gStyle->SetOptFit(1);
  	gStyle->SetOptTitle(0);
  	gStyle->SetTitleSize(0.05);

	if(dims == 2){
		// Setup histo
    	TH2F *h = new TH2F("h","3+3 #Chi^2;U_{e 4};U_{#mu 4}",1000,0.01,100.,1000,.01,100.);
		TH2F *h_99 = new TH2F("h99","",1000,0.0001,1.,1000,.01,100.);
		TH2F *h_90 = new TH2F("h90","",1000,0.0001,1.,1000,.01,100.);
    	h->GetXaxis()->SetTitleOffset(1.1);
    	h->GetYaxis()->SetTitleOffset(.8);
    	h->GetXaxis()->SetTitleFont(62);
    	h->GetYaxis()->SetTitleFont(62);
    	h->GetYaxis()->CenterTitle();
    	h->GetXaxis()->CenterTitle();
    	h->GetXaxis()->SetTitleSize(0.04);
    	h->GetXaxis()->SetLabelSize(0.04);
    	h->GetXaxis()->SetLabelOffset(0.001);
    	h->GetYaxis()->SetTitleSize(0.045);
    	h->GetYaxis()->SetLabelSize(0.04);
    	h->SetStats(kFALSE);

		c1->SetLogy();
		c1->SetLogx();
		chi2_99->SetMarkerStyle(7);
		chi2_90->SetMarkerStyle(7);
		chi2_99->SetFillColor(62);
    	chi2_90->SetFillColor(92);
		chi2_99->SetMarkerColor(62);
		chi2_90->SetMarkerColor(92);

		// Make overlay of paper plots for diag
		TGraph *overlay; Double_t x95[2100]; Double_t y95[2100];
		if(diag == 1){
			std::cout << " Diagnostic plot" << std::endl;
			std::string det;
			if(dataset == "karmen") det = "KARMEN";
			if(dataset == "ccfr") det = "CCFR";
			if(dataset == "cdhs") det = "CDHS";
			if(dataset == "lsnd") det = "LSND";
			if(dataset == "nomad") det = "nomad";
			if(dataset == "xsec") det = "XSEC";
			if(dataset == "bugey") det = "BUGEY";
			if(dataset == "numi") det = "NUMI";
			if(dataset == "mbnu") det = "MBNU";
			if(dataset == "mbnubar") det = "MBNUBAR";
			if(dataset == "mbnudis") det = "MBNUDIS";
			if(dataset == "gal") det = "GAL";
			if(dataset == "minos") det = "MINOS";
			std::string overFile = det + ".csv";
			std::cout << "Overlay: " << overFile << std::endl;

			if(overFile != ".csv"){
				ifstream file;
	    		file.open(overFile);
	    		for(int i = 0; i < 2100; i++){
	        		file >> x95[i];
					file >> y95[i];
				}
	    		file.close();
			}

			overlay = new TGraph(2100,x95,y95);
			overlay->SetFillColor(kAzure+1);
			overlay->SetMarkerColor(kAzure+1);
			overlay->SetMarkerStyle(7);
		}

		TLegend *leg = new TLegend(0.65,0.15,0.95,0.4);
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(62);
		leg->SetTextSize(0.04);
		if(diag == 0){
			leg->AddEntry(chi2_99,"99\% CL","f");
			leg->AddEntry(chi2_90,"90\% CL","f");
		}

		if(steriles == 3){
			h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{51}");
			h->GetXaxis()->SetLimits(.01,100.);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_99->Draw("m5*m5:m4*m4","","same");
			chi2_90->Draw("m5*m5:m4*m4","","same");
			leg->SetHeader("(3+3) Global Fit");
			leg->Draw();
        	c1->Print((plotOutput + "/" + dataset + "_dm251xdm241" + suffix + ".png").c_str());

			h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{61}");
			h->GetXaxis()->SetLimits(.01,100.);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_99->Draw("(m6*m6):(m4*m4)","","same");
			chi2_90->Draw("(m6*m6):(m4*m4)","","same");
			leg->SetHeader("(3+3) Global Fit");
			leg->Draw();
        	c1->Print(("plots/" + dataset + "_dm261xdm241" + suffix + ".png").c_str());
		}

		if(steriles == 2){
			h->SetTitle("#chi^{2} for 3+2 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{51}");
			h->GetXaxis()->SetLimits(.01,100.);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_99->Draw("m5*m5:m4*m4","","same");
			chi2_90->Draw("m5*m5:m4*m4","","same");
			leg->SetHeader("(3+2) Global Fit");
			leg->Draw();
        	c1->Print((plotOutput + "/" + dataset + "_3plus2_dm251xdm241" + suffix + ".png").c_str());
		}

		if(steriles == 1){
			if(raster == 0){
				if(type==0)	h->SetTitle("#chi^{2} for 3+1 Sterile Fits;sin^{2}(2#Theta_{#mue});#Deltam^{2}_{41}");
				if(type==0)	h->GetXaxis()->SetLimits(.0001,.1);
				if(type>0)	h->GetXaxis()->SetLimits(.0001,1.);
				h->GetYaxis()->SetLimits(.01,100.);
				h->Draw();

				// Now, draw overlay of the old one
				if(diag == 1){
					overlay->Draw("psame");
					leg->Draw();
				}
				else{
					if(type==0){
        				chi2_99->Draw("m4**2:4*um4**2*ue4**2","","same");
        				chi2_90->Draw("m4**2:4*um4**2*ue4**2","","same");
					}
					else if(type == 1){
						chi2_99->Draw("m4*m4:4*um4*um4*(1-um4*um4)","","same");
        				chi2_90->Draw("m4*m4:4*um4*um4*(1-um4*um4)","","same");
					}
					else if(type == 2){
						chi2_99->Draw("m4*m4:4*ue4*ue4*(1-ue4*ue4)","","same");
        				chi2_90->Draw("m4*m4:4*ue4*ue4*(1-ue4*ue4)","","same");
					}
					leg->SetHeader("(3+1) Global Fit");
					leg->Draw();
				}
        		c1->Print((plotOutput + "/" + dataset + "_3plus1_dm241xsinsq2t" + suffix + ".png").c_str());
			}
		}
	}

	// S'more plots
	// We'd also like to do a quick little analysis of 1d distributions, so let's throw those in as well.
	if(dims == 1){

		TCanvas * c2 = new TCanvas();
		c2->SetLogy();
		h_chi2->Draw();
		TLine *linechi2 = new TLine(chi2min,0,chi2min,1000);	linechi2->SetLineColor(2); 	linechi2->SetLineStyle(3); linechi2->SetLineWidth(3); 	linechi2->Draw();
		c2->Print((plotOutput + "/onedee_chi2_" + to_string(steriles) + ".png").c_str());
		h_ue4->Draw();
		TLine *lineue4 = new TLine(ue4_min,0,ue4_min,1000);	lineue4->SetLineColor(2); 	lineue4->SetLineStyle(3); lineue4->SetLineWidth(3); 	lineue4->Draw();
		c2->Print((plotOutput + "/onedee_ue4_" + to_string(steriles) + ".png").c_str());
		h_um4->Draw();
		TLine *lineum4 = new TLine(um4_min,0,um4_min,1000);	lineum4->SetLineColor(2); 	lineum4->SetLineStyle(3); lineum4->SetLineWidth(3); 	lineum4->Draw();
		c2->Print((plotOutput + "/onedee_um4_" + to_string(steriles) + ".png").c_str());
		if(steriles > 1){
			h_ue5->Draw();
			TLine *lineue5 = new TLine(ue5_min,0,ue5_min,1000);	lineue5->SetLineColor(2); 	lineue5->SetLineStyle(3); lineue5->SetLineWidth(3); 	lineue5->Draw();
			c2->Print((plotOutput + "/onedee_ue5_" + to_string(steriles) + ".png").c_str());
			h_um5->Draw();
			TLine *lineum5 = new TLine(um5_min,0,um5_min,1000);	lineum5->SetLineColor(2); 	lineum5->SetLineStyle(3); lineum5->SetLineWidth(3); 	lineum5->Draw();
			c2->Print((plotOutput + "/onedee_um5_" + to_string(steriles) + ".png").c_str());
			h_phi45->Draw();
			TLine *linephi45 = new TLine(phi45_min,0,phi45_min,1000);	linephi45->SetLineColor(2); 	linephi45->SetLineStyle(3); linephi45->SetLineWidth(3); 	linephi45->Draw();
			c2->Print((plotOutput + "/onedee_phi45_" + to_string(steriles) + ".png").c_str());
		}
		if(steriles > 2){
			h_ue6->Draw();
			TLine *lineue6 = new TLine(ue6_min,0,ue6_min,1000);	lineue6->SetLineColor(2); 	lineue6->SetLineStyle(3); lineue6->SetLineWidth(3); 	lineue6->Draw();
			c2->Print((plotOutput + "/onedee_ue6_" + to_string(steriles) + ".png").c_str());
			h_um6->Draw();
			TLine *lineum6 = new TLine(um6_min,0,um6_min,1000);	lineum6->SetLineColor(2); 	lineum6->SetLineStyle(3); lineum6->SetLineWidth(3); 	lineum6->Draw();
			c2->Print((plotOutput + "/onedee_um6_" + to_string(steriles) + ".png").c_str());
			h_phi46->Draw();
			TLine *linephi46 = new TLine(phi46_min,0,phi46_min,1000);	linephi46->SetLineColor(2); 	linephi46->SetLineStyle(3); linephi46->SetLineWidth(3); 	linephi46->Draw();
			c2->Print((plotOutput + "/onedee_phi46_" + to_string(steriles) + ".png").c_str());
			h_phi56->Draw();
			TLine *linephi56 = new TLine(phi56_min,0,phi56_min,1000);	linephi56->SetLineColor(2); 	linephi56->SetLineStyle(3); linephi56->SetLineWidth(3); 	linephi56->Draw();
			c2->Print((plotOutput + "/onedee_phi56_" + to_string(steriles) + ".png").c_str());
		}
		c2->SetLogx();
		h_m4->Draw();
		TLine *linem4 = new TLine(m4_min,0,m4_min,1000);	linem4->SetLineColor(2); 	linem4->SetLineStyle(3); linem4->SetLineWidth(3); 	linem4->Draw();
		c2->Print((plotOutput + "/onedee_m4_" + to_string(steriles) + ".png").c_str());
		if(steriles  > 1){
			h_m5->Draw();
			TLine *linem5 = new TLine(m5_min,0,m5_min,1000);	linem5->SetLineColor(2); 	linem5->SetLineStyle(3); linem5->SetLineWidth(3); 	linem5->Draw();
			c2->Print((plotOutput + "/onedee_m5_" + to_string(steriles) + ".png").c_str());
		}
		if(steriles > 2){
			h_m6->Draw();
			TLine *linem6 = new TLine(m6_min,0,m6_min,1000);	linem6->SetLineColor(2); 	linem6->SetLineStyle(3); linem6->SetLineWidth(3); 	linem6->Draw();
			c2->Print((plotOutput + "/onedee_m6_" + to_string(steriles) + ".png").c_str());
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
    globFit_plotter();
    return 0;
}
# endif

void plotter(){
    globFit_plotter();
    return;
}
