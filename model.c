#include "model.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TRandom.h"
#include <iostream>
#include "TLine.h"
#include "TH2D.h"
#include <vector>
#include "prob.h"

#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGFrame.h"
#include "TVirtualPad.h"


/*********************************************
* Data Struct: Same as Davio's output.
* -------------------------------------------
* ToDo: 
* Should be packed in a nTuple to me later, make sure remains consistant even if i change things. SO FAR DONE 
*
* search this for BUG 1 and find out why the hell outputting to canvas is broke and dependant on order of histograms?!
*
* get the oscillation probabilty right for e and f  DONE
*
* put an internal variable for e or mu into detector DONE
* *******************************************/




std::vector<double > SBN_spectrum::add_SBN_spectrum(SBN_spectrum other)
{
	std::vector<double > ans;
	std::vector<double > vsignal = SBN_spectrum::get_ninevector();
	std::vector<double > vbkg = other.get_ninevector();


	//First 3*N_e_bins are the SBND, muboone and icaraus e appearance bins so they are added to the background
	for(int k=0; k<3*(N_e_bins);k++)
	{
		ans.push_back(vbkg[k]+vsignal[k]);

	//	std::cout<<vbkg[k]<<" "<<vsignal[k]<<std::endl;
	}


	//Last 3*N_m_bins are the SBND, muboone and icaraus mu appearance bins so they are just the oscillated away ones.
	for(int k=3*(N_e_bins); k < 3*(N_e_bins+N_m_bins);k++)
	{

		ans.push_back(vsignal[k]);

	//	std::cout<<vbkg[k]<<" "<<vsignal[k]<<std::endl;

	}
	return ans;

}


SBN_spectrum::SBN_spectrum(struct neutrinoModel numodel){
	nullModel.zero();
	workingModel=numodel;	


	//Bit overboard in initilising these to be fair	
	sbnd_e.resize(N_e_bins);
	uboone_e.resize(N_e_bins);
	icarus_e.resize(N_e_bins);
	fill (sbnd_e.begin(),sbnd_e.end(),0.0); 
	fill (uboone_e.begin(),uboone_e.end(),0.0); 
	fill (icarus_e.begin(),icarus_e.end(),0.0); 

	sbnd_f.resize(N_e_bins);
	uboone_f.resize(N_e_bins);
	icarus_f.resize(N_e_bins);
	fill (sbnd_f.begin(),sbnd_f.end(),0.0); 
	fill (uboone_f.begin(),uboone_f.end(),0.0); 
	fill (icarus_f.begin(),icarus_f.end(),0.0); 

	sbnd_m.resize(N_m_bins);
	uboone_m.resize(N_m_bins);
	icarus_m.resize(N_m_bins);
	fill (sbnd_m.begin(),sbnd_m.end(),0.0); 
	fill (uboone_m.begin(),uboone_m.end(),0.0); 
	fill (icarus_m.begin(),icarus_m.end(),0.0); 




};

SBN_spectrum::SBN_spectrum(){
	nullModel.zero();
	workingModel=nullModel;
	
	sbnd_e.resize(N_e_bins);
	uboone_e.resize(N_e_bins);
	icarus_e.resize(N_e_bins);
	fill (sbnd_e.begin(),sbnd_e.end(),0.0); 
	fill (uboone_e.begin(),uboone_e.end(),0.0); 
	fill (icarus_e.begin(),icarus_e.end(),0.0); 

	sbnd_m.resize(N_m_bins);
	uboone_m.resize(N_m_bins);
	icarus_m.resize(N_m_bins);
	fill (sbnd_m.begin(),sbnd_m.end(),0.0); 
	fill (uboone_m.begin(),uboone_m.end(),0.0); 
	fill (icarus_m.begin(),icarus_m.end(),0.0);

/*	sbnd_m_bkg.resize(N_m_bins);
	uboone_m_bkg.resize(N_m_bins);
	icarus_m_bkg.resize(N_m_bins);
	fill (sbnd_m_bkg.begin(),sbnd_m_bkg.end(),0.0); 
	fill (uboone_m_bkg.begin(),uboone_m_bkg.end(),0.0); 
	fill (icarus_m_bkg.begin(),icarus_m_bkg.end(),0.0); */


	TFile *fSBND_app = new TFile("bkg_data/SBND_bkg_app.root");
	TFile *fSBND_dis= new TFile("bkg_data/SBND_bkg_dis.root");
	TFile *fuBooNE_app = new TFile("bkg_data/uBooNE_bkg_app.root");
	TFile *fuBooNE_dis= new TFile("bkg_data/uBooNE_bkg_dis.root");
	TFile *fICARUS_app = new TFile("bkg_data/ICARUS_bkg_app.root");
	TFile *fICARUS_dis= new TFile("bkg_data/ICARUS_bkg_dis.root");

	TH1D  SBND_app_intrinsic= *((TH1D*)fSBND_app->Get("SBND_intrinsic_app"));
	TH1D  SBND_app_muon= *((TH1D*)fSBND_app->Get("SBND_muon_app"));
	TH1D  SBND_app_photon=*( (TH1D*)fSBND_app->Get("SBND_photon_app"));
	TH1D  SBND_dis = *((TH1D*)fSBND_dis->Get("SBND_reco"));

	TH1D  uBooNE_app_intrinsic= *((TH1D*)fuBooNE_app->Get("UBOONE_intrinsic_app"));
	TH1D  uBooNE_app_muon= *((TH1D*)fuBooNE_app->Get("UBOONE_muon_app"));
	TH1D  uBooNE_app_photon=*( (TH1D*)fuBooNE_app->Get("UBOONE_photon_app"));
	TH1D  uBooNE_dis = *((TH1D*)fuBooNE_dis->Get("UBOONE_reco"));

	TH1D  ICARUS_app_intrinsic= *((TH1D*)fICARUS_app->Get("ICARUS_intrinsic_app"));
	TH1D  ICARUS_app_muon= *((TH1D*)fICARUS_app->Get("ICARUS_muon_app"));
	TH1D  ICARUS_app_photon=*( (TH1D*)fICARUS_app->Get("ICARUS_photon_app"));
	TH1D  ICARUS_dis = *((TH1D*)fICARUS_dis->Get("ICARUS_reco"));

	for(int i =0; i<N_e_bins; i++)
	{
		sbnd_e[i]=  SBND_app_intrinsic.GetBinContent(i+1)+SBND_app_muon.GetBinContent(i+1)+SBND_app_photon.GetBinContent(i+1);
		uboone_e[i]=  uBooNE_app_intrinsic.GetBinContent(i+1)+uBooNE_app_muon.GetBinContent(i+1)+uBooNE_app_photon.GetBinContent(i+1);
		icarus_e[i]=  ICARUS_app_intrinsic.GetBinContent(i+1)+ICARUS_app_muon.GetBinContent(i+1)+ICARUS_app_photon.GetBinContent(i+1);
	}

	for(int i=0; i < N_m_bins; i++)
	{
		sbnd_m[i]=  SBND_dis.GetBinContent(i+1);	
		uboone_m[i]=  uBooNE_dis.GetBinContent(i+1);	
		icarus_m[i]=  ICARUS_dis.GetBinContent(i+1);	
	}


}

void SBN_spectrum::vec_print(){
	std::cout<<"# bin_low  SBND   uBooNE   ICARUS "<<std::endl;
	std::cout<<"#******************** Oscillated nu_e ****************"<<std::endl;
	
	for(int k=0; k<sbnd_f.size(); k++)
	{
		std::cout<<e_bins[k]<<" "<<sbnd_f[k]<<" "<<uboone_f[k]<<" "<<icarus_f[k]<<std::endl;
	}


	std::cout<<"#******************** intrinsic nu_e ****************"<<std::endl;

	for(int k=0; k<sbnd_e.size(); k++)
	{
		std::cout<<e_bins[k]<<" "<<sbnd_e[k]<<" "<<uboone_e[k]<<" "<<icarus_e[k]<<std::endl;
	}

	std::cout<<"#******************** disap nu_mu ****************"<<std::endl;

	for(int k=0; k<sbnd_m.size(); k++)
	{
		std::cout<<mu_bins[k]<<" "<<sbnd_m[k]<<" "<<uboone_m[k]<<" "<<icarus_m[k]<<std::endl;
	}

}


const int SBN_spectrum::N_e_bins;
const int SBN_spectrum::N_m_bins;

constexpr double SBN_spectrum::mu_bins[N_m_bins+1];
constexpr double SBN_spectrum::e_bins[N_e_bins+1];



std::vector<double > SBN_spectrum::get_ninevector(){
			std::vector<double> ans = sbnd_f;
			ans.insert(std::end(ans), std::begin(uboone_f), std::end(uboone_f));
			ans.insert(std::end(ans), std::begin(icarus_f), std::end(icarus_f));
			ans.insert(std::end(ans), std::begin(sbnd_e), std::end(sbnd_e));
			ans.insert(std::end(ans), std::begin(uboone_e), std::end(uboone_e));
			ans.insert(std::end(ans), std::begin(icarus_e), std::end(icarus_e));
			ans.insert(std::end(ans), std::begin(sbnd_m), std::end(sbnd_m));
			ans.insert(std::end(ans), std::begin(uboone_m), std::end(uboone_m));
			ans.insert(std::end(ans), std::begin(icarus_m), std::end(icarus_m));
	
			return ans; 
};


std::vector<double > SBN_spectrum::get_sixvector(){

			std::vector<double > sbnd_tmp = sbnd_f;
			std::vector<double > uboone_tmp = uboone_f;
			std::vector<double > icarus_tmp = icarus_f;
			for(int k =0; k< sbnd_f.size(); k++)
			{

				sbnd_tmp[k] += sbnd_e[k];	
				uboone_tmp[k] += uboone_e[k];	
				icarus_tmp[k] += icarus_e[k];	
			}
	

			std::vector<double> ans = sbnd_tmp;
			ans.insert(std::end(ans), std::begin(uboone_tmp), std::end(uboone_tmp));
			ans.insert(std::end(ans), std::begin(icarus_tmp), std::end(icarus_tmp));
			ans.insert(std::end(ans), std::begin(sbnd_m), std::end(sbnd_m));
			ans.insert(std::end(ans), std::begin(uboone_m), std::end(uboone_m));
			ans.insert(std::end(ans), std::begin(icarus_m), std::end(icarus_m));
	
			return ans; 
};



void SBN_spectrum::update_model(struct neutrinoModel numodel){
	workingModel=numodel;
};

void SBN_spectrum::oscillate(){
//	uboone_e = BKG+ProbOscillate(MUBOONE_FLAG,E_FLAG);	
	//Will load the FULL_OSC, nu_E and nu_mu and eight them appropiately for each 
	//detector and each true L/E_v and so forth

	SBN_detector * SBND = new SBN_detector(DET_SBND);
	SBN_detector * UBOONE = new SBN_detector(DET_UBOONE);
	SBN_detector * ICARUS = new SBN_detector(DET_ICARUS);

	SBN_spectrum::fill_app(SBND);
	SBN_spectrum::fill_app(UBOONE);
	SBN_spectrum::fill_app(ICARUS);

	SBN_spectrum::fill_intrin(SBND);	
	SBN_spectrum::fill_intrin(UBOONE);	
	SBN_spectrum::fill_intrin(ICARUS);	

	//internal fiducial volume changes for muon run so lets take that into account
	bool mu_mode = true;
	SBN_detector * SBND_mu = new SBN_detector(DET_SBND, mu_mode);
	SBN_detector * UBOONE_mu = new SBN_detector(DET_UBOONE,mu_mode);
	SBN_detector * ICARUS_mu = new SBN_detector(DET_ICARUS,mu_mode);

	SBN_spectrum::fill_dis(SBND_mu);	
	SBN_spectrum::fill_dis(UBOONE_mu);	
	SBN_spectrum::fill_dis(ICARUS_mu);	


};

void SBN_spectrum::oscillate_sample(){
//	uboone_e = BKG+ProbOscillate(MUBOONE_FLAG,E_FLAG);	
	//Will load the FULL_OSC, nu_E and nu_mu and eight them appropiately for each 
	//detector and each true L/E_v and so forth

	SBN_detector * SBND = new SBN_detector(DET_SBND);
	SBN_detector * UBOONE = new SBN_detector(DET_UBOONE);
	SBN_detector * ICARUS = new SBN_detector(DET_ICARUS);

	SBN_spectrum::fill_app_sample(SBND);
	SBN_spectrum::fill_app_sample(UBOONE);
	SBN_spectrum::fill_app_sample(ICARUS);

	SBN_spectrum::fill_intrin_sample(SBND);	
	SBN_spectrum::fill_intrin_sample(UBOONE);	
	SBN_spectrum::fill_intrin_sample(ICARUS);	

	//internal fiducial volume changes for muon run so lets take that into account
	bool mu_mode = true;
	SBN_detector * SBND_mu = new SBN_detector(DET_SBND, mu_mode);
	SBN_detector * UBOONE_mu = new SBN_detector(DET_UBOONE,mu_mode);
	SBN_detector * ICARUS_mu = new SBN_detector(DET_ICARUS,mu_mode);

	SBN_spectrum::fill_dis_sample(SBND_mu);	
	SBN_spectrum::fill_dis_sample(UBOONE_mu);	
	SBN_spectrum::fill_dis_sample(ICARUS_mu);	


};


// REDUNDANT
void SBN_spectrum::fill_hists(){
	for(int i =0; i<uboone_e.size(); i++)
	{
		THuboone_e->SetBinContent(i+1,uboone_e[i]);
		THsbnd_e->SetBinContent(i+1,sbnd_e[i]);
		THicarus_e->SetBinContent(i+1,icarus_e[i]);
	}

	for(int i =0; i<uboone_m.size(); i++)
	{
		THuboone_m->SetBinContent(i+1,uboone_m[i]);
		THsbnd_m->SetBinContent(i+1,sbnd_m[i]);
		THicarus_m->SetBinContent(i+1,icarus_m[i]);
	}
};


// REDUNDANT
void SBN_spectrum::fill_vectors(){
	for(int i =0; i<THuboone_e->GetSize()-1; i++)
	{
		uboone_e.push_back( THuboone_e->GetBinContent(i+1));
		sbnd_e.push_back(    THsbnd_e->GetBinContent(i+1));
		icarus_e.push_back(  THicarus_e->GetBinContent(i+1));
	}

	for(int i=0; i<THuboone_m->GetSize()-1; i++)
	{
		uboone_m.push_back( THuboone_m->GetBinContent(i+1));
		sbnd_m.push_back(    THsbnd_m->GetBinContent(i+1));
		icarus_m.push_back(  THicarus_m->GetBinContent(i+1));
	}
};

	

int SBN_spectrum::fill_dis(SBN_detector * detector )
{


//	double DMSQ = workingModel.dm41Sq;
//	double S2TH = 1.0-4*pow(workingModel.Um[0],2)*(1- pow(workingModel.Um[0],2));


	TFile *fnudetector = new TFile(detector->fname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");

	TH1D H1nue("Muon Disapearance","Backgroun",N_m_bins,mu_bins);


	TRandom *rangen  = new TRandom();//initialize random generator	
 	rangen->SetSeed(8392.31);	
	
	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;

	double POTscaling      = detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume; 


	double Eff_em = 0.8*POTscaling*fiducialscaling;
	double Eff_ph = 0.06*POTscaling*fiducialscaling*0.192;
	double Eff_pi = 0.06*POTscaling*fiducialscaling*(1-0.9137);

	int pol = 0;
	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	int Ntest=0;

	double posX,posY,posZ;

	int Ncontained1 = 0;
	int Ncontained2 = 0;


	double Ep[20];
	double Epip[5];
	double Epim[5];
	double pdgo[12];
	double Eo[5];
	double Eph[5];
	double Epi0dph[10];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);

	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);
	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pl",pl);

	double vertex_pos[3] = {0,0,0};

	for(int i=0; i < tnudetector->GetEntries(); i++)
	{
	
	if(i%1000000==0){std::cout<<"Dis-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}



		tnudetector->GetEntry(i);
		
		detector->random_pos(rangen,vertex_pos); 


		//Is there a visible vertex and how much energy is there!
			
		
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
				//if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 

		//Check if we actually have a "visibe vertex"
		if(Ehad >= vertex_thresh)
		{
			vis_vertex = true;
		}
		
						
		/******************************************************************	
		 * CC Intrinsic Nu_mu
		 * ****************************************************************/		
				
			
		if((PDGnu==14||PDGnu==-14) && CC==1 )//&& Nph==0 && Npi0dph==0)
		{
			El_smear = smear_energy(El_true, MUsmear, rangen);

			if(true ) //El_smear >= EM_thresh)
			{
				Nsignal++;

				Enu_reco = Ehad + El_smear;


				double mag = sqrt(pl[0]*pl[0]+pl[1]*pl[1]+pl[2]*pl[2]);
				std::vector<double > dir;
			       	dir.push_back(pl[0]/mag);
			        dir.push_back(pl[1]/mag);
			        dir.push_back(pl[2]/mag);

			//	std::cout<<El_smear<<" "<<pl[0]<<" "<<pl[1]<<" "<<pl[2]<<std::endl;
			//
				double Lmu = muon_track_length(El_smear); 
				double il = 0;
				double observable_L = 0;
				
				double endpos[3] = {0,0,0};
				get_endpoint(vertex_pos,Lmu, pl, endpos);

				//double prob =Pmm(detector->baseline+(vertex_pos[2]/1000.0),El_smear,DMSQ,S2TH);
				double prob = workingModel.oscProb(2,2,El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
				
				Ncontained1++;
				if(detector->is_fully_contained(vertex_pos, endpos)){
						observable_L = Lmu;
						Ncontained2++;	
						if(observable_L > 50.0)
						{
					
							H1nue.Fill(Enu_reco ,weight*Eff_em*prob);
						}


				} 
				else
				{
					observable_L = detector->track_length_escape(vertex_pos,endpos);
					if(observable_L > 1000)
					{
						
						H1nue.Fill(Enu_reco ,weight*Eff_em*prob);

					}
				}





				
			} //end 200Mev smeared cut	


		}//end nu_mu cc cut
/************************************************************************************************
 *				NC Pi_pm mimicing
 * **********************************************************************************************/		
	if( ( !(Npim == 1) != !(Npip ==1))  && CC==1 )//&& Nph==0 && Npi0dph==0)
		{

			double Etrue=0;
			if(Npim==1){

				Etrue = Epip[0];
			} else 
			{
				Etrue = Epim[0];
			}
		


			El_smear = smear_energy(Etrue, MUsmear, rangen);

			if(true ) //El_smear >= EM_thresh)
			{
				//THcc_El->Fill(El_smear,weight*Eff_em);
				//THcc_El_true->Fill(El_true,weight*Eff_em);
				Nsignal++;

				Enu_reco = Ehad;


				if(pion_track_length(Etrue)<50)
				{
					H1nue.Fill(Enu_reco ,weight*Eff_em);
				}
			}
	
		}//End NC pi + loop
	} //end event loop
	char namei[200];
	sprintf(namei, "bkg_data/%s_dis.root",detector->name);

	TFile f(namei,"RECREATE");
	H1nue.Write();	

/* BUG 1 BUG 1
std::cout<<"just before canvas creation"<<std::endl; 
TCanvas c1("c1");
std::cout<<"just after canvas creation, before draw"<<std::endl; 

H1nue.Draw();

std::cout<<"just afer draw before update"<<std::endl; 
c1.Update();
  
std::cout<<"just afer update before save"<<std::endl; 
c1.SaveAs("test.png");

std::cout<<"done"<<std::endl; 
*/


	for(int i =0; i < N_m_bins; i++)
	{
	//	std::cout<<"det: "<<detector->identifier<<" "<<H1nue_reco_intrinsic.GetBinContent(i+1)<<std::endl;
		switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_m[i]= H1nue.GetBinContent(i+1);
				break;
			case DET_UBOONE:
				uboone_m[i]= H1nue.GetBinContent(i+1);
				break;
			case DET_ICARUS:
				icarus_m[i]= H1nue.GetBinContent(i+1);
				break;
		}


	}

}


/****************************************************************************
 *			Appearance for all SBN
 * *************************************************************************/

int SBN_spectrum::fill_app(SBN_detector * detector )
{

 //	SBN_detector SBND(400,400,500,2*183.5,370,405,110);
			
	//std::cout<<muon_track_length(1.2)<<" "<<muon_track_length(0.12)<<std::endl;


//	double DMSQ = workingModel.dm41Sq;
//	double S2TH = 4*pow(workingModel.Ue[0]*workingModel.Um[0],2.0);

	TFile *fnudetector = new TFile(detector->foscname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");
	

	TRandom *rangen    = new TRandom();//initialize random generator	
	rangen->SetSeed(9424.1);	

	TH1D  H_osc_1("on1","",N_e_bins,e_bins);
	

	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;



	double POTscaling =detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume;



	double Eff_em = 0.8*POTscaling*fiducialscaling;

	int pol = 0;
	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	int Ntest=0;

	double posX,posY,posZ;



	double Ep[20];
	double Epip[5];
	double Epim[5];
	double pdgo[12];
	double Eo[5];
	double Eph[5];
	double Epi0dph[10];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);

	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);
	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pl",pl);

	double vertex_pos[3] = {0,0,0};

	for(int i=0; i< tnudetector->GetEntries(); i++)
	{
	

		tnudetector->GetEntry(i);
		
		detector->random_pos(rangen,vertex_pos); 


		if(i%1000000==0){std::cout<<"App-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}


		//Is there a visible vertex and how much energy is there!
			
		
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
				//if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 

		//Check if we actually have a "visibe vertex"
		if(Ehad >= vertex_thresh)
		{
			vis_vertex = true;
		}
	
			
/************************************************************************************************
 *				CC Nu_e after oscillation! Lets fill a few
 * **********************************************************************************************/		
		
			
		if((PDGnu==12||PDGnu==-12) && CC==1 )//&& Nph==0 && Npi0dph==0)
		{
		El_smear = smear_energy(El_true, EMsmear, rangen);
		double prob = workingModel.oscProb(2,1,El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
			
	//	H_nonosc->Fill(El_smear+Ehad, weight*Eff_em);	
		H_osc_1.Fill(El_smear+Ehad, weight*Eff_em*prob);

	//	H_osc_2->Fill(El_smear+Ehad, weight*Eff_em*Pmue(110+vertex_pos[2]/1000.0,El_smear,2.0,0.75));	
	//	H_osc_3->Fill(El_smear+Ehad, weight*Eff_em*Pmue(110+vertex_pos[2]/1000.0,El_smear,20.0,0.75));	
		
			
		//	std::cout<<H_osc_1.GetBinContent(2)<<" "<<El_smear+Ehad<<std::endl;
		} //end 200Mev smeared cut	



	} //end event loop


	char namei[200];
	sprintf(namei, "bkg_data/%s_app.root",detector->name);

	TFile f(namei,"RECREATE");
	H_osc_1.Write();

	for(int i =0; i < N_e_bins; i++)
	{
	switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_f[i]= H_osc_1.GetBinContent(i+1);
				break;
			case DET_UBOONE:
				uboone_f[i]= H_osc_1.GetBinContent(i+1);		
				break;
			case DET_ICARUS:
				icarus_f[i]= H_osc_1.GetBinContent(i+1);
				break;
		}

	}


return 1;


};




/****************************************************************************
 *			Intrinsic-ness for all SBN
 * *************************************************************************/

int SBN_spectrum::fill_intrin(SBN_detector * detector )
{


	TRandom *rangen    = new TRandom();//initialize random generator	
	rangen->SetSeed(3268.94);	

	TFile *fnudetector = new TFile(detector->fname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");

	TH1D THcc_El("detector_cc_el","",N_e_bins,e_bins);
	TH1D THcc_El_true("detector_cc_el_true","",N_e_bins,e_bins);
	TH1D THintrin_nue ("detector_nu_e","",N_e_bins,e_bins);
	TH1D H1nue_reco_intrinsic("detector_intrinsic_app","",N_e_bins,e_bins);
	TH1D H1nue_reco_muon("detector_muon_app","",N_e_bins,e_bins);
	TH1D H1nue_reco_photon("detector_photon_app","",N_e_bins,e_bins);

	TH1D H_check_muon("check muon","",25,0.2,0.65);
	TH1D H_check_muon2("check muon2","",25,0.10,0.65);


	TH1D H_muon_track("muon track length","",40,0,800);
	TH1D H_muon_track_contained("muon track lengthcont","",40,0,800);


	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;

	double POTscaling =detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume;


	double Eff_em = 0.8*POTscaling*fiducialscaling;
	double Eff_ph = 0.06*POTscaling*fiducialscaling;
	double Eff_pi = 0.06*POTscaling*fiducialscaling;
	double Eff_ver = 1.0;
	double Eff_in = 0.9137;
	double Eff_out = (1.0-Eff_in); 




	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	double Ntest=0.0;

	double vertex_pos[3] = {0,0,0};

	double Ep[100];
	double Epip[100];
	double Epim[100];
	int    pdgo[100];
	double Eo[100];
	double Eph[100];
	double Epi0dph[100];
	double pph[100][3];
	double ppi0dph[100][3];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);
	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);

	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pph",pph);
	tnudetector->SetBranchAddress("ppi0dph",ppi0dph);
	tnudetector->SetBranchAddress("pl",pl);

	for(int i=0; i< tnudetector->GetEntries(); i++)
	{

		tnudetector->GetEntry(i);
	
	if(i%1000000==0){std::cout<<"Intrinsic-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}

		
		detector->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 
		
		//Is there a visible vertex and how much energy is there!
				
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
			//	if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 



		//Check if we actually have a "visibe vertex"
		if(Ehad > vertex_thresh)
		{
			vis_vertex = true;
		} else 
		{
			vis_vertex = false;
		}
	
			
/************************************************************************************************
 *				CC Intrinsic Nu_e 
 * **********************************************************************************************/		
		
			
		if((PDGnu==12||PDGnu==-12) && CC== 1 && Nph==0 && Npi0dph==0)
		{
			El_smear = smear_energy(El_true, EMsmear, rangen);

			if(El_smear >= EM_thresh)
			{


		double prob = workingModel.oscProb(1,1,El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));

				THcc_El.Fill(El_smear, weight*Eff_em);
				THcc_El_true.Fill(El_true, weight*Eff_em);
				Nsignal++;

				Enu_reco = El_smear + Ehad;
				H1nue_reco_intrinsic.Fill(Enu_reco, weight*Eff_em*prob);
				
			} //end 200Mev smeared cut	


		}//end nu_e cc cut
/************************************************************************************************
 *				NC photon bkg 
 * **********************************************************************************************/	
		if(NC== 1 && (Nph!=0 || Npi0dph!=0)) //BEGIN NC 1 gamma part
		{

			if(Nph == 1 || Npi0dph == 1) //single photon NOT from pion
			{
				double Eph_smeared = 0;
				double Lph =0; 

				if(Nph == 1){	
					Eph_smeared = smear_energy(Eph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Eph[0],rangen); 

				} else {
					Eph_smeared = smear_energy(Epi0dph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Epi0dph[0],rangen); 
				}


				if(Eph_smeared >= EM_thresh)
				{

					if((vis_vertex && Lph < 3.0) || !vis_vertex)
					{
					  	H1nue_reco_photon.Fill(Eph_smeared + Ehad ,weight*Eff_ph);
					}
				
				}




			}//end single photon cut

			if(Npi0dph >= 3){Ntest=Ntest+weight*Eff_ph;}

			if(Nph == 2 || Npi0dph == 2) //double pion photons
			{
				double E1 =0;	
				double E2 =0;	
				double L1 = 0.0; 
				double L2 = 0.0;	
				double p1norm = 0;
				double p2norm = 0;
				double conv1[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};
				double conv2[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};


				if(Nph == 2){	
					 E1 = smear_energy(Eph[0],EMsmear,rangen);
					 E2 = smear_energy(Eph[1],EMsmear,rangen);
					 L1 =  photon_conversion_length(Eph[0],rangen);
			        	 L2 =  photon_conversion_length(Eph[1],rangen);
					 p1norm = sqrt(pow(pph[0][0],2)+pow(pph[0][1],2)+pow(pph[0][2],2));
					 p2norm = sqrt(pow(pph[1][0],2)+pow(pph[1][1],2)+pow(pph[1][2],2));

					 conv1[0] += L1*pph[0][0]/p1norm;
					 conv1[1] += L1*pph[0][1]/p1norm;
					 conv1[2] += L1*pph[0][2]/p1norm;
					 
					 conv2[0] += L1*pph[1][0]/p2norm;
					 conv2[1] += L1*pph[1][1]/p2norm;
					 conv2[2] += L1*pph[1][2]/p2norm;
				
				} else {
					E1 = smear_energy(Epi0dph[0],EMsmear,rangen);
					E2 = smear_energy(Epi0dph[1],EMsmear,rangen);
					L1 =  photon_conversion_length(Epi0dph[0],rangen);
			        	L2 =  photon_conversion_length(Epi0dph[1],rangen);
				       	p1norm = sqrt(pow(ppi0dph[0][0],2)+pow(ppi0dph[0][1],2)+pow(ppi0dph[0][2],2));
				        p2norm = sqrt(pow(ppi0dph[1][0],2)+pow(ppi0dph[1][1],2)+pow(ppi0dph[1][2],2));
				
					 conv1[0] += L1*ppi0dph[0][0]/p1norm;
					 conv1[1] += L1*ppi0dph[0][1]/p1norm;
					 conv1[2] += L1*ppi0dph[0][2]/p1norm;
					 
					 conv2[0] += L1*ppi0dph[1][0]/p2norm;
					 conv2[1] += L1*ppi0dph[1][1]/p2norm;
					 conv2[2] += L1*ppi0dph[1][2]/p2norm;
				}

			
				if( detector->is_fiducial(conv1) && !detector->is_fiducial(conv2))
				{
					// Then we only have a single photon, photon 1

					if(E1 >= EM_thresh)
					{
					     	if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_photon.Fill(E1 + Ehad ,weight*Eff_ph);
						}
					}
				} 
				else if( detector->is_fiducial(conv2) && !detector->is_fiducial(conv1))
				{
					// Then we only have a single photon, photon 2
					if(E2 >= EM_thresh)
					{
	
					     	if(( vis_vertex && L2 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_photon.Fill(E2+ Ehad ,weight*Eff_ph);
						}
					}

				} 
				else if( detector->is_active(conv2) && detector->is_active(conv1))
				{

		
					if(E1 >= EM_thresh && E2 < 0.1 && detector->is_fiducial(conv1)) 
					{
					     		if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
							{	
						     		H1nue_reco_photon.Fill(E1 + Ehad  , weight*Eff_ph);
							}
 	
					}
					else if(E2 >= EM_thresh && E1 < 0.1 && detector->is_fiducial(conv2)) 
					{
						     	if((vis_vertex && L2 < 3.0) || !vis_vertex)
							{	
					     			H1nue_reco_photon.Fill(E2 + Ehad , weight*Eff_ph);
							}
					}
				}
		
				


			}//end single pion 

		
		}//End photon NC 


/************************************************************************************************
 *				CC Muon Nu_mu + Gamma 
 * **********************************************************************************************/	
		if(( PDGnu==14 || PDGnu==-14) && CC==1)
		{
		double fudge = 0.0;


		double El_smear = smear_energy(El_true,MUsmear,rangen);
		double observable_L = 0;
		double track_L = muon_track_length(El_true);
		double endpoint[3] = {0,0,0}; //will store position of end of muon track (may be outside detector volume)


		get_endpoint(vertex_pos,track_L, pl, endpoint);


		if(detector->is_fully_contained(vertex_pos, endpoint)){

			observable_L = track_L;

			H_muon_track_contained.Fill(track_L,weight);
		} 
		else 
		{
			observable_L = detector->track_length_escape(vertex_pos,endpoint);
		}

		H_muon_track.Fill(observable_L,weight);




		if( observable_L < 100.0)	//It is track length, ones can be less if they are above vertex
		{

			H_check_muon2.Fill(El_true,weight);

			//First Check if there is a EM shower
	
			if(Nph == 1 || Npi0dph == 1) //single photon NOT from pion
			{
		
				double Eph_smeared = 0;
				double Lph =0; 
				if(Nph == 1){	
					Eph_smeared = smear_energy(Eph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Eph[0],rangen); 

				} else {
					Eph_smeared = smear_energy(Epi0dph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Epi0dph[0],rangen); 
				}

				if(Eph_smeared >= EM_thresh)
				{

					if((vis_vertex && Lph < 3.0) || !vis_vertex)
					{
					  	H1nue_reco_muon.Fill(Eph_smeared + Ehad +El_smear -fudge ,weight*Eff_ph);
						H_check_muon.Fill(Eph_smeared + Ehad +El_smear -fudge ,weight*Eff_ph);

					}
				
				}


			}//end single photon cut


			if(Nph == 2 || Npi0dph == 2) //double pion photons
			{
				double E1 =0;	
				double E2 =0;	
				double L1 = 0.0; 
				double L2 = 0.0;	
				double p1norm = 0;
				double p2norm = 0;
				double conv1[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};
				double conv2[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};

				if(Nph == 2){	
					 E1 = smear_energy(Eph[0],EMsmear,rangen);
					 E2 = smear_energy(Eph[1],EMsmear,rangen);
					 L1 =  photon_conversion_length(Eph[0],rangen);
			        	 L2 =  photon_conversion_length(Eph[1],rangen);
					 p1norm = sqrt(pow(pph[0][0],2)+pow(pph[0][1],2)+pow(pph[0][2],2));
					 p2norm = sqrt(pow(pph[1][0],2)+pow(pph[1][1],2)+pow(pph[1][2],2));

					 conv1[0] += L1*pph[0][0]/p1norm;
					 conv1[1] += L1*pph[0][1]/p1norm;
					 conv1[2] += L1*pph[0][2]/p1norm;
					 
					 conv2[0] += L1*pph[1][0]/p2norm;
					 conv2[1] += L1*pph[1][1]/p2norm;
					 conv2[2] += L1*pph[1][2]/p2norm;
				
				} else {
					E1 = smear_energy(Epi0dph[0],EMsmear,rangen);
					E2 = smear_energy(Epi0dph[1],EMsmear,rangen);
					L1 =  photon_conversion_length(Epi0dph[0],rangen);
			        	L2 =  photon_conversion_length(Epi0dph[1],rangen);
				       	p1norm = sqrt(pow(ppi0dph[0][0],2)+pow(ppi0dph[0][1],2)+pow(ppi0dph[0][2],2));
				        p2norm = sqrt(pow(ppi0dph[1][0],2)+pow(ppi0dph[1][1],2)+pow(ppi0dph[1][2],2));
				
					 conv1[0] += L1*ppi0dph[0][0]/p1norm;
					 conv1[1] += L1*ppi0dph[0][1]/p1norm;
					 conv1[2] += L1*ppi0dph[0][2]/p1norm;
					 
					 conv2[0] += L1*ppi0dph[1][0]/p2norm;
					 conv2[1] += L1*ppi0dph[1][1]/p2norm;
					 conv2[2] += L1*ppi0dph[1][2]/p2norm;
				}

			
				if( detector->is_fiducial(conv1) && !detector->is_fiducial(conv2))
				{
					// Then we only have a single photon, photon 1

					if(E1 >= EM_thresh)
					{
					     	if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_muon.Fill(E1 + Ehad+El_smear -fudge,weight*Eff_ph);
							H_check_muon.Fill(E1 + Ehad +El_smear -fudge ,weight*Eff_ph);
						}
					}
				} 
				else if( detector->is_fiducial(conv2) && !detector->is_fiducial(conv1))
				{
					// Then we only have a single photon, photon 2
					if(E2 >= EM_thresh)
					{
	
					     	if(( vis_vertex && L2 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_muon.Fill(E2+ Ehad + El_smear -fudge ,weight*Eff_ph);
							H_check_muon.Fill(E2 + Ehad +El_smear -fudge ,weight*Eff_ph);
						}
					}

				} 
				else if( detector->is_active(conv2) && detector->is_active(conv1))
				{

		
					if(E1 >= EM_thresh && E2 < 0.1 && detector->is_fiducial(conv1)) 
					{
					     		if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
							{	
						     		H1nue_reco_muon.Fill(E1 + Ehad +El_smear -fudge  , weight*Eff_ph);
								H_check_muon.Fill(E1 + Ehad +El_smear -fudge ,weight*Eff_ph);
							}
 	
					}
					else if(E2 >= EM_thresh && E1 < 0.1 && detector->is_fiducial(conv2)) 
					{
						     	if((vis_vertex && L2 < 3.0) || !vis_vertex)
							{	
					     			H1nue_reco_muon.Fill(E2 + Ehad +El_smear- fudge, weight*Eff_ph);
								H_check_muon.Fill(E2 + Ehad +El_smear -fudge ,weight*Eff_ph);
							}
					}
				}
		
				


				


			}//end single pion 





		}}//End Muon Cut





	}



	char namei[200];
	sprintf(namei, "bkg_data/%s_intrin.root",detector->name);

	TFile f(namei,"RECREATE");

	H1nue_reco_intrinsic.Write();
	H1nue_reco_muon.Write();
	H1nue_reco_photon.Write();

	
	
	for(int i =0; i < N_e_bins; i++)
	{
	switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);
				break;
			case DET_UBOONE:
			        uboone_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);				    
			       	break;
			case DET_ICARUS:
				icarus_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);
				break;
		}

	}


};



//**********************************************************************************************
//
//
//
//
//
//**********************************************************************************************

	

int SBN_spectrum::fill_dis_sample(SBN_detector * detector )
{


//	double DMSQ = workingModel.dm41Sq;
//	double S2TH = 1.0-4*pow(workingModel.Um[0],2)*(1- pow(workingModel.Um[0],2));


	TFile *fnudetector = new TFile(detector->fname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");


	TH1D H1nue_sin("dis_muon_sin","Backgroun",N_m_bins,mu_bins);
	TH1D H1nue_sinsq("dis_muon_sinsq","Backgroun",N_m_bins,mu_bins);

	TRandom *rangen  = new TRandom();//initialize random generator	
 	rangen->SetSeed(8392.31);	
	
	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;

	double POTscaling      = detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume; 


	double Eff_em = 0.8*POTscaling*fiducialscaling;
	double Eff_ph = 0.06*POTscaling*fiducialscaling*0.192;
	double Eff_pi = 0.06*POTscaling*fiducialscaling*(1-0.9137);

	int pol = 0;
	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	int Ntest=0;

	double posX,posY,posZ;

	int Ncontained1 = 0;
	int Ncontained2 = 0;


	double Ep[20];
	double Epip[5];
	double Epim[5];
	double pdgo[12];
	double Eo[5];
	double Eph[5];
	double Epi0dph[10];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);

	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);
	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pl",pl);

	double vertex_pos[3] = {0,0,0};

	for(int i=0; i < tnudetector->GetEntries(); i++)
	{
	
	if(i%1000000==0){std::cout<<"Dis-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}



		tnudetector->GetEntry(i);
		
		detector->random_pos(rangen,vertex_pos); 


		//Is there a visible vertex and how much energy is there!
			
		
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
				//if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 

		//Check if we actually have a "visibe vertex"
		if(Ehad >= vertex_thresh)
		{
			vis_vertex = true;
		}
		
						
		/******************************************************************	
		 * CC Intrinsic Nu_mu
		 * ****************************************************************/		
				
			
		if((PDGnu==14||PDGnu==-14) && CC==1 )//&& Nph==0 && Npi0dph==0)
		{
			El_smear = smear_energy(El_true, MUsmear, rangen);

			if(true ) //El_smear >= EM_thresh)
			{
				Nsignal++;

				Enu_reco = Ehad + El_smear;


				double mag = sqrt(pl[0]*pl[0]+pl[1]*pl[1]+pl[2]*pl[2]);
				std::vector<double > dir;
			       	dir.push_back(pl[0]/mag);
			        dir.push_back(pl[1]/mag);
			        dir.push_back(pl[2]/mag);

			//	std::cout<<El_smear<<" "<<pl[0]<<" "<<pl[1]<<" "<<pl[2]<<std::endl;
			//
				double Lmu = muon_track_length(El_smear); 
				double il = 0;
				double observable_L = 0;
				
				double endpos[3] = {0,0,0};
				get_endpoint(vertex_pos,Lmu, pl, endpos);

				//double prob =Pmm(detector->baseline+(vertex_pos[2]/1000.0),El_smear,DMSQ,S2TH);
				
				double prob = workingModel.oscProbSin(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
				double probsq = workingModel.oscProbSinSq(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
				Ncontained1++;
				if(detector->is_fully_contained(vertex_pos, endpos)){
						observable_L = Lmu;
						Ncontained2++;	
						if(observable_L > 50.0)
						{
					
							H1nue_sin.Fill(Enu_reco ,weight*Eff_em*prob);
							H1nue_sinsq.Fill(Enu_reco ,weight*Eff_em*probsq);
						}


				} 
				else
				{
					observable_L = detector->track_length_escape(vertex_pos,endpos);
					if(observable_L > 1000)
					{
						
						H1nue_sin.Fill(Enu_reco ,weight*Eff_em*prob);
						H1nue_sinsq.Fill(Enu_reco ,weight*Eff_em*probsq);

					}
				}





				
			} //end 200Mev smeared cut	


		}//end nu_mu cc cut
/************************************************************************************************
 *				NC Pi_pm mimicing
 * **********************************************************************************************/		
	if( ( !(Npim == 1) != !(Npip ==1))  && CC==1 )//&& Nph==0 && Npi0dph==0)
		{

			double Etrue=0;
			if(Npim==1){

				Etrue = Epip[0];
			} else 
			{
				Etrue = Epim[0];
			}
		


			El_smear = smear_energy(Etrue, MUsmear, rangen);

			if(true ) //El_smear >= EM_thresh)
			{
				//THcc_El->Fill(El_smear,weight*Eff_em);
				//THcc_El_true->Fill(El_true,weight*Eff_em);
				Nsignal++;

				Enu_reco = Ehad;


				if(pion_track_length(Etrue)<50)
				{
					H1nue_sin.Fill(Enu_reco ,weight*Eff_em);
					H1nue_sinsq.Fill(Enu_reco ,weight*Eff_em);
				}
			}
	
		}//End NC pi + loop
	} //end event loop
	char namei[200];
	sprintf(namei, "bkg_data/%s_%.2f.root",detector->name,workingModel.dm41Sq);

	TFile f(namei,"UPDATE");

	H1nue_sin.Write();	
	H1nue_sinsq.Write();

/* BUG 1 BUG 1
std::cout<<"just before canvas creation"<<std::endl; 
TCanvas c1("c1");
std::cout<<"just after canvas creation, before draw"<<std::endl; 

H1nue.Draw();

std::cout<<"just afer draw before update"<<std::endl; 
c1.Update();
  
std::cout<<"just afer update before save"<<std::endl; 
c1.SaveAs("test.png");

std::cout<<"done"<<std::endl; 
*/

/*
	for(int i =0; i < N_m_bins; i++)
	{
	//	std::cout<<"det: "<<detector->identifier<<" "<<H1nue_reco_intrinsic.GetBinContent(i+1)<<std::endl;
		switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_m[i]= H1nue.GetBinContent(i+1);
				break;
			case DET_UBOONE:
				uboone_m[i]= H1nue.GetBinContent(i+1);
				break;
			case DET_ICARUS:
				icarus_m[i]= H1nue.GetBinContent(i+1);
				break;
		}


	}

	*/

}


/****************************************************************************
 *			Appearance for all SBN
 * *************************************************************************/

int SBN_spectrum::fill_app_sample(SBN_detector * detector )
{

 //	SBN_detector SBND(400,400,500,2*183.5,370,405,110);
			
	//std::cout<<muon_track_length(1.2)<<" "<<muon_track_length(0.12)<<std::endl;


//	double DMSQ = workingModel.dm41Sq;
//	double S2TH = 4*pow(workingModel.Ue[0]*workingModel.Um[0],2.0);

	TFile *fnudetector = new TFile(detector->foscname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");
	

	TRandom *rangen    = new TRandom();//initialize random generator	
	rangen->SetSeed(9424.1);	

	
	TH1D  H_osc_1_sin("fullosc_nue_sin","",N_e_bins,e_bins);
	TH1D  H_osc_1_sinsq("fullosc_nue_sinsq","",N_e_bins,e_bins);

	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;



	double POTscaling =detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume;



	double Eff_em = 0.8*POTscaling*fiducialscaling;

	int pol = 0;
	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	int Ntest=0;

	double posX,posY,posZ;



	double Ep[20];
	double Epip[5];
	double Epim[5];
	double pdgo[12];
	double Eo[5];
	double Eph[5];
	double Epi0dph[10];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);

	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);
	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pl",pl);

	double vertex_pos[3] = {0,0,0};

	for(int i=0; i< tnudetector->GetEntries(); i++)
	{
	

		tnudetector->GetEntry(i);
		
		detector->random_pos(rangen,vertex_pos); 


		if(i%1000000==0){std::cout<<"App-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}


		//Is there a visible vertex and how much energy is there!
			
		
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
				//if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 

		//Check if we actually have a "visibe vertex"
		if(Ehad >= vertex_thresh)
		{
			vis_vertex = true;
		}
	
			
/************************************************************************************************
 *				CC Nu_e after oscillation! Lets fill a few
 * **********************************************************************************************/		
		
			
		if((PDGnu==12||PDGnu==-12) && CC==1 )//&& Nph==0 && Npi0dph==0)
		{
		El_smear = smear_energy(El_true, EMsmear, rangen);
			
		double prob = workingModel.oscProbSin(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
		double probsq = workingModel.oscProbSinSq(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
	//	H_nonosc->Fill(El_smear+Ehad, weight*Eff_em);	

		H_osc_1_sin.Fill(El_smear+Ehad, weight*Eff_em*prob);
		H_osc_1_sinsq.Fill(El_smear+Ehad, weight*Eff_em*probsq);
	//	H_osc_2->Fill(El_smear+Ehad, weight*Eff_em*Pmue(110+vertex_pos[2]/1000.0,El_smear,2.0,0.75));	
	//	H_osc_3->Fill(El_smear+Ehad, weight*Eff_em*Pmue(110+vertex_pos[2]/1000.0,El_smear,20.0,0.75));	
		
			
		//	std::cout<<H_osc_1.GetBinContent(2)<<" "<<El_smear+Ehad<<std::endl;
		} //end 200Mev smeared cut	



	} //end event loop


	char namei[200];
	sprintf(namei, "bkg_data/%s_%.2f.root",detector->name,workingModel.dm41Sq);

	TFile f(namei,"UPDATE");
	H_osc_1_sin.Write();
	H_osc_1_sinsq.Write();
	
	/*
	for(int i =0; i < N_e_bins; i++)
	{
	switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_f[i]= H_osc_1.GetBinContent(i+1);
				break;
			case DET_UBOONE:
				uboone_f[i]= H_osc_1.GetBinContent(i+1);		
				break;
			case DET_ICARUS:
				icarus_f[i]= H_osc_1.GetBinContent(i+1);
				break;
		}

	}
*/

return 1;


};




/****************************************************************************
 *			Intrinsic-ness for all SBN
 * *************************************************************************/

int SBN_spectrum::fill_intrin_sample(SBN_detector * detector )
{


	TRandom *rangen    = new TRandom();//initialize random generator	
	rangen->SetSeed(3268.94);	

	TFile *fnudetector = new TFile(detector->fname);
	TTree *tnudetector = (TTree*)fnudetector->Get("mainTree");

	TH1D THcc_El("detector_cc_el","",N_e_bins,e_bins);
	TH1D THcc_El_true("detector_cc_el_true","",N_e_bins,e_bins);
	TH1D THintrin_nue ("detector_nu_e","",N_e_bins,e_bins);
	
	TH1D H1nue_reco_intrinsic_sin("nue_intrin_sin","",N_e_bins,e_bins);
	TH1D H1nue_reco_muon_sin("nue_muon_sin","",N_e_bins,e_bins);
	TH1D H1nue_reco_photon_sin("nue_photon_sin","",N_e_bins,e_bins);
	TH1D H1nue_reco_intrinsic_sinsq("nue_intrin_sinsq","",N_e_bins,e_bins);
	TH1D H1nue_reco_muon_sinsq("nue_muon_sinsq","",N_e_bins,e_bins);
	TH1D H1nue_reco_photon_sinsq("nue_photon_sinsq","",N_e_bins,e_bins);



	TH1D H_check_muon("check muon","",25,0.2,0.65);
	TH1D H_check_muon2("check muon2","",25,0.10,0.65);


	TH1D H_muon_track("muon track length","",40,0,800);
	TH1D H_muon_track_contained("muon track lengthcont","",40,0,800);


	int Nnue = 0;
	int Nnuebar=0;
	int Nnumu = 0;
	int Nnumubar=0;
	int Nsignal = 0;

	int NpOBS = 0;
	int NpipOBS = 0;
	int NpimOBS = 0;

	double POTscaling =detector->potmodifier*66.0; 
	double fiducialscaling = detector->f_volume/detector->volume;


	double Eff_em = 0.8*POTscaling*fiducialscaling;
	double Eff_ph = 0.06*POTscaling*fiducialscaling;
	double Eff_pi = 0.06*POTscaling*fiducialscaling;
	double Eff_ver = 1.0;
	double Eff_in = 0.9137;
	double Eff_out = (1.0-Eff_in); 




	double Enu;
	double El_true;
	double El_smear;
	int PDGnu;
	double weight;
	int Np;
	int Npip;
	int Npim;
	int CC;
	int NC;
	int Nph;
	int Npi0dph;
	int No;

	double Ntest=0.0;

	double vertex_pos[3] = {0,0,0};

	double Ep[100];
	double Epip[100];
	double Epim[100];
	int    pdgo[100];
	double Eo[100];
	double Eph[100];
	double Epi0dph[100];
	double pph[100][3];
	double ppi0dph[100][3];
	double pl[3];

	tnudetector->SetBranchAddress("Enu",&Enu);
	tnudetector->SetBranchAddress("El",&El_true);
	tnudetector->SetBranchAddress("PDGnu",&PDGnu);
	tnudetector->SetBranchAddress("weight",&weight);
	tnudetector->SetBranchAddress("Np",&Np);
	tnudetector->SetBranchAddress("Npip",&Npip);
	tnudetector->SetBranchAddress("Npim",&Npim);
	tnudetector->SetBranchAddress("Nph",&Nph);
	tnudetector->SetBranchAddress("No",&No);
	tnudetector->SetBranchAddress("Npi0dph",&Npi0dph);
	tnudetector->SetBranchAddress("CC",&CC);
	tnudetector->SetBranchAddress("NC",&NC);

	tnudetector->SetBranchAddress("Ep",Ep);
	tnudetector->SetBranchAddress("Eph",Eph);
	tnudetector->SetBranchAddress("Epim",Epim);
	tnudetector->SetBranchAddress("Epip",Epip);
	tnudetector->SetBranchAddress("Epi0dph",Epi0dph);
	tnudetector->SetBranchAddress("Eo",Eo);
	tnudetector->SetBranchAddress("pdgo",pdgo);
	tnudetector->SetBranchAddress("pph",pph);
	tnudetector->SetBranchAddress("ppi0dph",ppi0dph);
	tnudetector->SetBranchAddress("pl",pl);

	for(int i=0; i< tnudetector->GetEntries(); i++)
	{

		tnudetector->GetEntry(i);
	
	if(i%1000000==0){std::cout<<"Intrinsic-Det: "<<detector->identifier<<" #: "<<i<<std::endl;}

		
		detector->random_pos(rangen,vertex_pos); //Assign a random position for vertex in detector 
		
		//Is there a visible vertex and how much energy is there!
				
		double Enu_reco = 0;
		double Ehad=0;
		bool vis_vertex = false;


		if(Np!=0){
			double p_kin_true = 0;
			double p_kin_smeared = 0;

			for(int j=0; j<Np; j++)
			{
				p_kin_true = Ep[j]-MPROTON;
				p_kin_smeared = smear_energy(p_kin_true,psmear,rangen);
				if(p_kin_smeared>p_thresh)
				{
					Ehad += p_kin_smeared;
				}
			}
		} //end proton addition 

		if(Npip!=0){
			double pip_kin_true = 0;
			double pip_kin_smeared = 0;
			for(int j=0; j<Npip;j++)
			{
				pip_kin_true = Epip[j]-MPION;
				pip_kin_smeared = smear_energy(pip_kin_true,pismear,rangen);
				if(pip_kin_smeared>pip_thresh)
				{
					Ehad += pip_kin_smeared+MPION;
				}
						

			} 
		}//end pi+addition

		if(Npim!=0){
			double pim_kin_true = 0;
			double pim_kin_smeared = 0;
			for(int j=0; j<Npim;j++)
			{
				pim_kin_true = Epim[j]-MPION;
				pim_kin_smeared = smear_energy(pim_kin_true,pismear,rangen);
				if(pim_kin_smeared > pim_thresh)
				{
					Ehad += pim_kin_smeared+MPION;
				}
				
			} 
		}//end piminus addition

		if(No!=0){
					
			for(int j=0; j<No;j++)
			{
				Ehad += smear_energy(Eo[j],pismear,rangen);
			
			//	if(pdgo[j]!=0){	std::cout<<pdgo[j]<<std::endl;	}
			} 
		}//end Other 



		//Check if we actually have a "visibe vertex"
		if(Ehad > vertex_thresh)
		{
			vis_vertex = true;
		} else 
		{
			vis_vertex = false;
		}
	
			
/************************************************************************************************
 *				CC Intrinsic Nu_e 
 * **********************************************************************************************/		
		
			
		if((PDGnu==12||PDGnu==-12) && CC== 1 && Nph==0 && Npi0dph==0)
		{
			El_smear = smear_energy(El_true, EMsmear, rangen);

			if(El_smear >= EM_thresh)
			{


			double prob = workingModel.oscProbSin(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));
			double probsq = workingModel.oscProbSinSq(El_smear,0.001*(detector->baseline+(vertex_pos[2]/1000.0)));

				THcc_El.Fill(El_smear, weight*Eff_em);
				THcc_El_true.Fill(El_true, weight*Eff_em);
				Nsignal++;

				Enu_reco = El_smear + Ehad;
				H1nue_reco_intrinsic_sin.Fill(Enu_reco, weight*Eff_em*prob);
				H1nue_reco_intrinsic_sinsq.Fill(Enu_reco, weight*Eff_em*probsq);
				
			} //end 200Mev smeared cut	


		}//end nu_e cc cut
/************************************************************************************************
 *				NC photon bkg 
 * **********************************************************************************************/	
		if(NC== 1 && (Nph!=0 || Npi0dph!=0)) //BEGIN NC 1 gamma part
		{

			if(Nph == 1 || Npi0dph == 1) //single photon NOT from pion
			{
				double Eph_smeared = 0;
				double Lph =0; 

				if(Nph == 1){	
					Eph_smeared = smear_energy(Eph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Eph[0],rangen); 

				} else {
					Eph_smeared = smear_energy(Epi0dph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Epi0dph[0],rangen); 
				}


				if(Eph_smeared >= EM_thresh)
				{

					if((vis_vertex && Lph < 3.0) || !vis_vertex)
					{
					  	H1nue_reco_photon_sin.Fill(Eph_smeared + Ehad ,weight*Eff_ph);
					  	H1nue_reco_photon_sinsq.Fill(Eph_smeared + Ehad ,weight*Eff_ph);
					}
				
				}




			}//end single photon cut

			if(Npi0dph >= 3){Ntest=Ntest+weight*Eff_ph;}

			if(Nph == 2 || Npi0dph == 2) //double pion photons
			{
				double E1 =0;	
				double E2 =0;	
				double L1 = 0.0; 
				double L2 = 0.0;	
				double p1norm = 0;
				double p2norm = 0;
				double conv1[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};
				double conv2[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};


				if(Nph == 2){	
					 E1 = smear_energy(Eph[0],EMsmear,rangen);
					 E2 = smear_energy(Eph[1],EMsmear,rangen);
					 L1 =  photon_conversion_length(Eph[0],rangen);
			        	 L2 =  photon_conversion_length(Eph[1],rangen);
					 p1norm = sqrt(pow(pph[0][0],2)+pow(pph[0][1],2)+pow(pph[0][2],2));
					 p2norm = sqrt(pow(pph[1][0],2)+pow(pph[1][1],2)+pow(pph[1][2],2));

					 conv1[0] += L1*pph[0][0]/p1norm;
					 conv1[1] += L1*pph[0][1]/p1norm;
					 conv1[2] += L1*pph[0][2]/p1norm;
					 
					 conv2[0] += L1*pph[1][0]/p2norm;
					 conv2[1] += L1*pph[1][1]/p2norm;
					 conv2[2] += L1*pph[1][2]/p2norm;
				
				} else {
					E1 = smear_energy(Epi0dph[0],EMsmear,rangen);
					E2 = smear_energy(Epi0dph[1],EMsmear,rangen);
					L1 =  photon_conversion_length(Epi0dph[0],rangen);
			        	L2 =  photon_conversion_length(Epi0dph[1],rangen);
				       	p1norm = sqrt(pow(ppi0dph[0][0],2)+pow(ppi0dph[0][1],2)+pow(ppi0dph[0][2],2));
				        p2norm = sqrt(pow(ppi0dph[1][0],2)+pow(ppi0dph[1][1],2)+pow(ppi0dph[1][2],2));
				
					 conv1[0] += L1*ppi0dph[0][0]/p1norm;
					 conv1[1] += L1*ppi0dph[0][1]/p1norm;
					 conv1[2] += L1*ppi0dph[0][2]/p1norm;
					 
					 conv2[0] += L1*ppi0dph[1][0]/p2norm;
					 conv2[1] += L1*ppi0dph[1][1]/p2norm;
					 conv2[2] += L1*ppi0dph[1][2]/p2norm;
				}

			
				if( detector->is_fiducial(conv1) && !detector->is_fiducial(conv2))
				{
					// Then we only have a single photon, photon 1

					if(E1 >= EM_thresh)
					{
					     	if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_photon_sin.Fill(E1 + Ehad ,weight*Eff_ph);
						  	H1nue_reco_photon_sinsq.Fill(E1 + Ehad ,weight*Eff_ph);
						}
					}
				} 
				else if( detector->is_fiducial(conv2) && !detector->is_fiducial(conv1))
				{
					// Then we only have a single photon, photon 2
					if(E2 >= EM_thresh)
					{
	
					     	if(( vis_vertex && L2 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_photon_sin.Fill(E2+ Ehad ,weight*Eff_ph);

						  	H1nue_reco_photon_sinsq.Fill(E2+ Ehad ,weight*Eff_ph);
						}
					}

				} 
				else if( detector->is_active(conv2) && detector->is_active(conv1))
				{

		
					if(E1 >= EM_thresh && E2 < 0.1 && detector->is_fiducial(conv1)) 
					{
					     		if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
							{	
						     		H1nue_reco_photon_sin.Fill(E1 + Ehad  , weight*Eff_ph);
						     		H1nue_reco_photon_sinsq.Fill(E1 + Ehad  , weight*Eff_ph);
							}
 	
					}
					else if(E2 >= EM_thresh && E1 < 0.1 && detector->is_fiducial(conv2)) 
					{
						     	if((vis_vertex && L2 < 3.0) || !vis_vertex)
							{	

					     			H1nue_reco_photon_sin.Fill(E2 + Ehad , weight*Eff_ph);
					     			H1nue_reco_photon_sinsq.Fill(E2 + Ehad , weight*Eff_ph);
							}
					}
				}
		
				


			}//end single pion 

		
		}//End photon NC 


/************************************************************************************************
 *				CC Muon Nu_mu + Gamma 
 * **********************************************************************************************/	
		if(( PDGnu==14 || PDGnu==-14) && CC==1)
		{
		double fudge = 0.0;


		double El_smear = smear_energy(El_true,MUsmear,rangen);
		double observable_L = 0;
		double track_L = muon_track_length(El_true);
		double endpoint[3] = {0,0,0}; //will store position of end of muon track (may be outside detector volume)


		get_endpoint(vertex_pos,track_L, pl, endpoint);


		if(detector->is_fully_contained(vertex_pos, endpoint)){

			observable_L = track_L;

			H_muon_track_contained.Fill(track_L,weight);
		} 
		else 
		{
			observable_L = detector->track_length_escape(vertex_pos,endpoint);
		}

		H_muon_track.Fill(observable_L,weight);




		if( observable_L < 100.0)	//It is track length, ones can be less if they are above vertex
		{

			H_check_muon2.Fill(El_true,weight);

			//First Check if there is a EM shower
	
			if(Nph == 1 || Npi0dph == 1) //single photon NOT from pion
			{
		
				double Eph_smeared = 0;
				double Lph =0; 
				if(Nph == 1){	
					Eph_smeared = smear_energy(Eph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Eph[0],rangen); 

				} else {
					Eph_smeared = smear_energy(Epi0dph[0],EMsmear,rangen);
				        Lph = photon_conversion_length(Epi0dph[0],rangen); 
				}

				if(Eph_smeared >= EM_thresh)
				{

					if((vis_vertex && Lph < 3.0) || !vis_vertex)
					{
					  	H1nue_reco_muon_sin.Fill(Eph_smeared + Ehad +El_smear -fudge ,weight*Eff_ph);
					  	H1nue_reco_muon_sinsq.Fill(Eph_smeared + Ehad +El_smear -fudge ,weight*Eff_ph);
						
						H_check_muon.Fill(Eph_smeared + Ehad +El_smear -fudge ,weight*Eff_ph);

					}
				
				}


			}//end single photon cut


			if(Nph == 2 || Npi0dph == 2) //double pion photons
			{
				double E1 =0;	
				double E2 =0;	
				double L1 = 0.0; 
				double L2 = 0.0;	
				double p1norm = 0;
				double p2norm = 0;
				double conv1[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};
				double conv2[3] = {vertex_pos[0],vertex_pos[1],vertex_pos[2]};

				if(Nph == 2){	
					 E1 = smear_energy(Eph[0],EMsmear,rangen);
					 E2 = smear_energy(Eph[1],EMsmear,rangen);
					 L1 =  photon_conversion_length(Eph[0],rangen);
			        	 L2 =  photon_conversion_length(Eph[1],rangen);
					 p1norm = sqrt(pow(pph[0][0],2)+pow(pph[0][1],2)+pow(pph[0][2],2));
					 p2norm = sqrt(pow(pph[1][0],2)+pow(pph[1][1],2)+pow(pph[1][2],2));

					 conv1[0] += L1*pph[0][0]/p1norm;
					 conv1[1] += L1*pph[0][1]/p1norm;
					 conv1[2] += L1*pph[0][2]/p1norm;
					 
					 conv2[0] += L1*pph[1][0]/p2norm;
					 conv2[1] += L1*pph[1][1]/p2norm;
					 conv2[2] += L1*pph[1][2]/p2norm;
				
				} else {
					E1 = smear_energy(Epi0dph[0],EMsmear,rangen);
					E2 = smear_energy(Epi0dph[1],EMsmear,rangen);
					L1 =  photon_conversion_length(Epi0dph[0],rangen);
			        	L2 =  photon_conversion_length(Epi0dph[1],rangen);
				       	p1norm = sqrt(pow(ppi0dph[0][0],2)+pow(ppi0dph[0][1],2)+pow(ppi0dph[0][2],2));
				        p2norm = sqrt(pow(ppi0dph[1][0],2)+pow(ppi0dph[1][1],2)+pow(ppi0dph[1][2],2));
				
					 conv1[0] += L1*ppi0dph[0][0]/p1norm;
					 conv1[1] += L1*ppi0dph[0][1]/p1norm;
					 conv1[2] += L1*ppi0dph[0][2]/p1norm;
					 
					 conv2[0] += L1*ppi0dph[1][0]/p2norm;
					 conv2[1] += L1*ppi0dph[1][1]/p2norm;
					 conv2[2] += L1*ppi0dph[1][2]/p2norm;
				}

			
				if( detector->is_fiducial(conv1) && !detector->is_fiducial(conv2))
				{
					// Then we only have a single photon, photon 1

					if(E1 >= EM_thresh)
					{
					     	if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_muon_sin.Fill(E1 + Ehad+El_smear -fudge,weight*Eff_ph);
						  	H1nue_reco_muon_sinsq.Fill(E1 + Ehad+El_smear -fudge,weight*Eff_ph);
							H_check_muon.Fill(E1 + Ehad +El_smear -fudge ,weight*Eff_ph);
						}
					}
				} 
				else if( detector->is_fiducial(conv2) && !detector->is_fiducial(conv1))
				{
					// Then we only have a single photon, photon 2
					if(E2 >= EM_thresh)
					{
	
					     	if(( vis_vertex && L2 < 3.0 ) || !vis_vertex)
						{
						  	H1nue_reco_muon_sin.Fill(E2+ Ehad + El_smear -fudge ,weight*Eff_ph);
						  	H1nue_reco_muon_sinsq.Fill(E2+ Ehad + El_smear -fudge ,weight*Eff_ph);
							H_check_muon.Fill(E2 + Ehad +El_smear -fudge ,weight*Eff_ph);
						}
					}

				} 
				else if( detector->is_active(conv2) && detector->is_active(conv1))
				{

		
					if(E1 >= EM_thresh && E2 < 0.1 && detector->is_fiducial(conv1)) 
					{
					     		if(( vis_vertex && L1 < 3.0 ) || !vis_vertex)
							{	
						     		H1nue_reco_muon_sin.Fill(E1 + Ehad +El_smear -fudge  , weight*Eff_ph);
						     		H1nue_reco_muon_sinsq.Fill(E1 + Ehad +El_smear -fudge  , weight*Eff_ph);
								H_check_muon.Fill(E1 + Ehad +El_smear -fudge ,weight*Eff_ph);
							}
 	
					}
					else if(E2 >= EM_thresh && E1 < 0.1 && detector->is_fiducial(conv2)) 
					{
						     	if((vis_vertex && L2 < 3.0) || !vis_vertex)
							{	
					     			H1nue_reco_muon_sin.Fill(E2 + Ehad +El_smear- fudge, weight*Eff_ph);
					     			H1nue_reco_muon_sinsq.Fill(E2 + Ehad +El_smear- fudge, weight*Eff_ph);
								H_check_muon.Fill(E2 + Ehad +El_smear -fudge ,weight*Eff_ph);
							}
					}
				}
		
				


				


			}//end single pion 





		}}//End Muon Cut





	}



	char namei[200];
	sprintf(namei, "bkg_data/%s_%.2f.root",detector->name,workingModel.dm41Sq);

	TFile f(namei,"UPDATE");

	H1nue_reco_intrinsic_sin.Write();
	H1nue_reco_muon_sin.Write();
	H1nue_reco_photon_sin.Write();
	H1nue_reco_intrinsic_sinsq.Write();
	H1nue_reco_muon_sinsq.Write();
	H1nue_reco_photon_sinsq.Write();




/*	
	for(int i =0; i < N_e_bins; i++)
	{
	switch (detector->identifier)
		{
			case DET_SBND:
				sbnd_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);
				break;
			case DET_UBOONE:
			        uboone_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);				    
			       	break;
			case DET_ICARUS:
				icarus_e[i]= H1nue_reco_intrinsic.GetBinContent(i+1)+H1nue_reco_muon.GetBinContent(i+1)+H1nue_reco_photon.GetBinContent(i+1);
				break;
		}

	}
*/


};



