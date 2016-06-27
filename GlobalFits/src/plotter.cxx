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

int globFit_plotter(){

	procOptLoc = "/Users/dcianci/Physics/SBN_3plusN/GlobalFits/inputs/";
    procOpt();

	std::cout << "TYPE: " << type << std::endl;

    std::cout << "Loading ntuple files..." << std::endl;
	std::string infile;
	if(discretized == 0)	infile = output + Form("/nt_3%i_",steriles) + dataset + ".root";
	if(discretized == 1)	infile = output + Form("/nt_3%i_",steriles) + dataset + "_processed.root";
	TString inputFile = infile;
	std::cout << "Infile: " << infile << std::endl;
	TFile *f = new TFile(inputFile);
	TNtuple *chi2_99;
	TNtuple *chi2_90;
	TNtuple *chi2_95;

	if(discretized == 0){
		chi2_99 = (TNtuple*)(f->Get("chi2_99"));
		chi2_90 = (TNtuple*)(f->Get("chi2_90"));
		chi2_95 = (TNtuple*)(f->Get("chi2_95"));
		suffix = "";
	}
	if(discretized == 1){
		chi2_99 = (TNtuple*)(f->Get("chi2_99_pr"));
		chi2_90 = (TNtuple*)(f->Get("chi2_90_pr"));
		suffix = "_disc";
	}

	// Find chi2Min
	chi2_90->SetBranchAddress("chi2",&chi2);
	chi2_90->SetBranchAddress("m4",&m4);
	chi2_90->SetBranchAddress("ue4",&ue4);
	chi2_90->SetBranchAddress("um4",&um4);
	chi2_90->SetBranchAddress("m5",&m5);
	chi2_90->SetBranchAddress("ue5",&ue5);
	chi2_90->SetBranchAddress("um5",&um5);
	chi2_90->SetBranchAddress("m6",&m6);
	chi2_90->SetBranchAddress("ue6",&ue6);
	chi2_90->SetBranchAddress("um6",&um6);
	chi2_90->SetBranchAddress("phi45",&phi45);
	chi2_90->SetBranchAddress("phi46",&phi46);
	chi2_90->SetBranchAddress("phi56",&phi56);
	float chi2min = 3000.f;
    for(int i = 0; i < chi2_90->GetEntries(); i++){
        chi2_90->GetEntry(i);
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

	TCanvas *c1 = new TCanvas("c1");

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
		if(diag == 1){
			chi2_95->SetMarkerStyle(7);
			chi2_95->SetFillColor(kRed+3);
			chi2_95->SetMarkerColor(kRed+3);
		}

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

		TLegend *leg = new TLegend(0.1,0.1,0.35,0.3);
		leg->SetFillStyle(0);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(62);
		leg->SetTextSize(0.03);
		if(diag == 0){
			leg->AddEntry(chi2_99,"99%% CL","f");
			leg->AddEntry(chi2_90,"90%% CL","f");
		}
		if(diag == 1){
			leg->AddEntry(chi2_95,"Davio's fits","f");
			leg->AddEntry(overlay,"Paper fits","f");
		}


		if(steriles == 3){
			h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{51}");
			h->GetXaxis()->SetLimits(.01,100.);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_99->Draw("m5*m5:m4*m4","","same");
			chi2_90->Draw("m5*m5:m4*m4","","same");
			leg->Draw();
        	c1->Print((plotOutput + "/" + dataset + "_dm251xdm241" + suffix + ".png").c_str());

			h->SetTitle("#chi^{2} for 3+3 Sterile Fits;#Delta m^{2}_{41};#Delta m^{2}_{61}");
			h->GetXaxis()->SetLimits(.01,100.);
			h->GetYaxis()->SetLimits(.01,100.);
			h->Draw();

			chi2_99->Draw("(m6*m6):(m4*m4)","","same");
			chi2_90->Draw("(m6*m6):(m4*m4)","","same");
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
			leg->Draw();
        	c1->Print((plotOutput + "/" + dataset + "_3plus2_dm251xdm241" + suffix + ".png").c_str());
		}

		if(steriles == 1){
			if(raster == 0){
				if(type==0)	h->SetTitle("#chi^{2} for 3+1 Sterile Fits;sin^{2}(2#Theta_{e#mu});#Delta m^{2}_{41}");
				if(type==0)	h->GetXaxis()->SetLimits(.0001,.1);
				if(type>0)	h->GetXaxis()->SetLimits(.0001,1.);
				h->GetYaxis()->SetLimits(.01,100.);
				h->Draw();

				// Now, draw overlay of the old one
				if(diag == 1){
					chi2_95->Draw("dm2:sin22th","","same");
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
					leg->Draw();
				}
        		c1->Print((plotOutput + "/" + dataset + "_3plus1_dm241xsinsq2t" + suffix + ".png").c_str());
			}
			if(raster == 1){
				if(type==0)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{e#mu});#Delta m^{2}_{41}");
				if(type==1)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{#mu#mu});#Delta m^{2}_{41}");
				if(type==2)	h->SetTitle("95%%CL for 3+1 Sterile Fits;sin^{2}(2#Theta_{ee});#Delta m^{2}_{41}");
				if(type==0)	h->GetXaxis()->SetLimits(.0001,.1);
				if(type>0)	h->GetXaxis()->SetLimits(.0001,1.);
				h->GetYaxis()->SetLimits(.01,100.);
				h->Draw();

				if(diag == 1)
					overlay->Draw("psame");
				chi2_95->Draw("dm2:sin22th","","same");
				leg->Draw();
				c1->Print((plotOutput + "/" + dataset + "_3plus1_dm241xsinsq2t_raster" + suffix + ".png").c_str());
			}
		}
	}

	if(dims == 3 && steriles == 1){

		// Shit, so the plots we want are dm2 vs sin22thmm, vs sin22thme and dm2 vs ue4, vs um4
		TH3F *h = new TH3F("h1","3+1 #Chi^{2};U_{#mu 4};U_{e 4};#Delta m^{2}_{41}",1000,0.0001,1.,1000,0.0001,1.,1000,0.01,100.);
		TH3F *h3d1_99 = new TH3F("h2","",1000,0.0001,1.,1000,0.0001,1.,1000,0.01,100.);
		TH3F *h3d1_90 = new TH3F("h3","",1000,0.0001,1.,1000,0.0001,1.,1000,0.01,100.);
		TH3F *h3d2_99 = new TH3F("h4 ","",1000,0.,.5,1000,0.,.5,1000,0.01,100.);
		TH3F *h3d2_90 = new TH3F("h5","",1000,0.,.5,1000,0.,.5,1000,0.01,100.);

		h->GetXaxis()->SetTitleOffset(1.4);
		h->GetYaxis()->SetTitleOffset(1.6);
		h->GetZaxis()->SetTitleOffset(.9);
		h->GetXaxis()->SetTitleFont(62);
		h->GetYaxis()->SetTitleFont(62);
		h->GetZaxis()->SetTitleFont(62);
		h->GetYaxis()->CenterTitle();
		h->GetXaxis()->CenterTitle();
		h->GetZaxis()->CenterTitle();
		h->GetXaxis()->SetTitleSize(0.04);
		h->GetXaxis()->SetLabelSize(0.035);
		h->GetYaxis()->SetTitleSize(0.04);
		h->GetYaxis()->SetLabelSize(0.035);
		h->GetYaxis()->SetTitleSize(0.04);
		h->GetYaxis()->SetLabelSize(0.035);
    	h->SetStats(kFALSE);

		c1->SetLogz();
		h3d1_99->SetMarkerStyle(7);			h3d2_99->SetMarkerStyle(7);
		h3d1_90->SetMarkerStyle(7);			h3d2_90->SetMarkerStyle(7);
		h3d1_99->SetFillColor(kBlue);		h3d2_99->SetFillColor(kBlue);
    	h3d1_90->SetFillColor(kMagenta);	h3d2_90->SetFillColor(kMagenta);
		h3d1_99->SetMarkerColor(kBlue);		h3d2_99->SetMarkerColor(kBlue);
		h3d1_90->SetMarkerColor(kMagenta);	h3d2_90->SetMarkerColor(kMagenta);

		// Now fill the histos:
		// First, 99%
		float m4, ue4, um4;
		chi2_99->SetBranchAddress("m4",&m4);
		chi2_99->SetBranchAddress("ue4",&ue4);
		chi2_99->SetBranchAddress("um4",&um4);
		float mstep = TMath::Log10(100./.01)/float(1000);
		float ustep = (.5)/float(1000);
		float sinstep = TMath::Log10(1./.0001)/float(1000);
		float sinsmm, sinsme;
		for(int i = 0; i < chi2_99->GetEntries(); i++){
			chi2_99->GetEntry(i);
			sinsme = 4*pow(um4,2)*pow(ue4,2); 	sinsmm = 4*pow(um4,2)*(1-pow(um4,2));
			h3d1_99->SetBinContent(ceil(TMath::Log10(sinsme/.0001)/sinstep),ceil(TMath::Log10(sinsmm/.0001)/sinstep),ceil(TMath::Log10(pow(m4,2)/.01)/mstep),1);
			h3d2_99->SetBinContent(ceil(um4/ustep),ceil(ue4/ustep),ceil(TMath::Log10(pow(m4,2)/.01)/mstep),1);
		}
		// now 90
		chi2_90->SetBranchAddress("m4",&m4);
		chi2_90->SetBranchAddress("ue4",&ue4);
		chi2_90->SetBranchAddress("um4",&um4);
		for(int i = 0; i < chi2_90->GetEntries(); i++){
			chi2_90->GetEntry(i);
			sinsme = 4*pow(um4,2)*pow(ue4,2); 	sinsmm = 4*pow(um4,2)*(1-pow(um4,2));
			h3d1_90->SetBinContent(ceil(TMath::Log10(sinsme/.0001)/sinstep),ceil(TMath::Log10(sinsmm/.0001)/sinstep),ceil(TMath::Log10(pow(m4,2)/.01)/mstep),1);
			h3d2_90->SetBinContent(ceil(um4/ustep),ceil(ue4/ustep),ceil(TMath::Log10(pow(m4,2)/.01)/mstep),1);
		}

		//h->Draw();
		h3d2_99->Draw("p");
		h3d2_90->Draw("psame");
		c1->Print((plotOutput + "/" + dataset + "_threedee_ues" + suffix + ".png").c_str());

		h->SetTitle("3+1 #Chi^{2};sin^{2}(2#Theta_{e#mu});sin^{2}(2#Theta_{#mu#mu});#Delta m^{2}_{41}");
		c1->SetLogx();
		c1->SetLogy();
		h->Draw();
		h3d1_99->Draw("psame");
		h3d1_90->Draw("psame");
		c1->Print((plotOutput + "/" + dataset + "_threedee_sinsqs" + suffix + ".png").c_str());
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
