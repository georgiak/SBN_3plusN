#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "model.h"
#include "correlation.h"
#include "TRandom.h"
#include "prob.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


/*************************************************************
 *************************************************************
 *Hardcoded Defines -> poor practice, will wrap intp SBN_spectrum soon
 actually no here is probably best.
 ************************************************************
 ************************************************************/

#define N_m_bins 19
#define N_e_bins 11
#define N_dets 3
#define N_e_spectra 7
#define N_m_spectra 2

#define no_argument 0
#define required_argument 1
#define optional_argument 2

/*************************************************************
 *************************************************************
 *		BEGIN Main::sbnfit.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

gSystem->Load("libTree");


// Just a few flags to control program flow.
bool fit_flag = false;
bool verbose_flag = false;
bool test_flag=false;
bool bkg_flag= false;
bool sample_flag = false;
bool cov_flag = false;
bool sens_flag=false;
bool dis_flag = false;
bool app_flag = false;
bool comb_flag = false;
bool both_flag = false;
bool unit_flag = false;
bool fraction_flag = false;

int  sens_num=1;
double dm41 = -1.0;

double in_dm = 0;
double in_ue4 = 0;
double in_um4=0;

bool stat_only = false;
int dis_which = 1;
int num_ster = 0;
int index; 
int iarg = 0;
opterr=1;
/*************************************************************
 *************************************************************
 *		Command Line Argument Reading
 ************************************************************
 ************************************************************/

const struct option longopts[] = 
{
	{"fit", 		no_argument, 		0, 'F'},
	{"bkg",			no_argument,		0, 'B'},
	{"test", 		no_argument, 		0, 'T'},
	{"mass",		required_argument,	0, 'm'},
	{"ue4", 		required_argument,	0, 'e'},
	{"um4"	,		required_argument,	0, 'u'},
	{"help",		no_argument, 		0, 'h'},
	{"verbose",		no_argument, 		0, 'v'},
	{"sensitivity",		required_argument,	0, 'S'},
	{"stat-only",		no_argument,		0, 'l'},
	{"sample",		no_argument,		0, 's'},
	{"cov",			no_argument, 		0, 'c'},
	{"dis",			required_argument, 	0, 'd'},
	{"unitary",		no_argument,		0, 'n'},
	{"fraction",		required_argument,	0, 'f'},
	{"app",			no_argument,		0, 'a'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "d:alf:nue:m:svhS:cFTB", longopts, &index);

	switch(iarg)
	{
		case 'F':
			fit_flag = true;
			//mS = strtof(optarg,NULL);
			break;
		case 'B':
			bkg_flag = true;
			break;
		case 'n':
			unit_flag = true;
			break;
		case 'T':
			test_flag = true;
			break;	
		case 'v':
			verbose_flag = true;
			break;
		case 's':
			sample_flag = true;
			break;
		case 'l':
			stat_only = true;
			break;
		case 'd':
			dis_flag = true;
			dis_which =strtof(optarg,NULL);
			break;
		case 'a':
			app_flag = true;
			break;
		case 'm':
			in_dm  = strtof(optarg,NULL);
			break;
		case 'f':
			fraction_flag = true;
			num_ster = strtof(optarg,NULL);
			break;
		case 'u':
			in_um4  = strtof(optarg,NULL);
			break;
		case 'e':
			in_ue4  = strtof(optarg,NULL);
			break;
		case 'c':
			cov_flag = true;
			break;
		case 'S':
			sens_flag = true;
			sens_num  = strtof(optarg,NULL);
			break;
		case '?':
		case 'h':
			std::cout<<"Allowed arguments:"<<std::endl;
			std::cout<<"\t-F\t--fit\t\tRun a single SBN fit. Used in conjuction with --mass, --ue4, --um4"<<std::endl;
			std::cout<<"\t-m\t--mass\t\tSet 3p1 dm41 value"<<std::endl;
			std::cout<<"\t-u\t--um4\t\tSet 3p1 Um4 value"<<std::endl;
			std::cout<<"\t-e\t--ue4\t\tSet 3p1 Ue4 value"<<std::endl;
			std::cout<<"\t-T\t--test\t\trun SBN test code"<<std::endl;
			std::cout<<"\t-S\t--sensitivity\t\trun a full sensitivity fit. Required argument, number of steriles. Run with -a -d"<<std::endl;	
			std::cout<<"\t-f\t--fraction\t\tRun fraction analysis, takes 1 argument for 3+N"<<std::endl;
			std::cout<<"\t-d\t--app\t\tRun app only sensitivity"<<std::endl;
			std::cout<<"\t-a\t--dis\t\tRun dis only sensitivity, takes 1 arg 0 for e-dis 1 for mudis"<<std::endl;
			std::cout<<"\t-l\t--stat-only\t\t Run with no systematics, ony stats"<<std::endl;
			std::cout<<"\t-B\t--bkg\t\trun bkg generating test code"<<std::endl;
			std::cout<<"\t-h\t--help\t\tDisplays this help message"<<std::endl;
			std::cout<<"\t-s\t--sample\t\t generate all the sin and sin^2 samples"<<std::endl;
			std::cout<<"\t-d\t\t\tRequired Argument. Creates a sin and sin^2 frequency ntuples for a dmsq."<<std::endl;
			std::cout<<"\t-v\t\t\tVerbose run, mostly debugging"<<std::endl;	
			std::cout<<"\t-n\t\t\tUnitary run"<<std::endl;
			return 0;
	}

}

if(app_flag&&dis_flag){both_flag = true; app_flag=false; dis_flag=false;}



if(unit_flag){

	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;


	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
	bkgspec.sbnd_e_dirt[1]=44*1/5;
	bkgspec.uboone_e_dirt[0]= 47*4/5;
	bkgspec.uboone_e_dirt[1]=47*1/5;
	bkgspec.icarus_e_dirt[0]= 67*4/5;
	bkgspec.icarus_e_dirt[1]=67*1/5;
	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;
	

	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();

		int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */
		TMatrixT <double> M(matrix_size,matrix_size);
		TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
		TMatrixT <double> McI(matrix_size_c, matrix_size_c);


		int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int contMsize = (N_e_bins+N_m_bins)*N_dets;

		TMatrixT <double> Msys(bigMsize,bigMsize);
		sys_fill(Msys,usedetsys);

		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back[i]*back[j];
			}
		}




		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back);

		//stat_only =false;
		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		//if(stat_only){
		//	Mtotal =  Mstat;
		//} else {
		//	Mtotal = Msys+Mstat;
		//}
		Mtotal = Mstat+Msys;

		TMatrixT<double > Mctotal(contMsize,contMsize);
		contract_signal2(Mtotal,Mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		McI = Mctotal.Invert(&invdet);


		for(double uuee = 1; uuee > 1-0.15; uuee-=0.01){
		for(double uuem = 0; uuem <0.05; uuem+=0.002){//0.005
       		for(double uumm = 1; uumm >1-0.08; uumm-=0.005){//0.02
				//double uuee =0;// rangen->Uniform(0,1);
				//double uumm =0;// rangen->Uniform(0,1);
			//	double uuee = 0.3; 
			//	double uumm = 0;
			//	double uumm = 1;
			//	double uuee = 1;
				//double uuem =0;// rangen->Uniform(0,1);
				double uume =uuem;// rangen->Uniform(0,1);

				neutrinoModel unitModel(0.0,0.0,0.0);
				unitModel.zero();
				unitModel.UUmm=uumm;
				unitModel.UUee=uuee;
				unitModel.UUem=uuem;
				unitModel.UUme=uume;

				SBN_spectrum UnitSpec(unitModel);
	
				UnitSpec.load_unit(ICARUS);
				UnitSpec.load_unit(SBND);
				UnitSpec.load_unit(UBOONE);

					
				UnitSpec.sbnd_e_dirt[1]=44*1/5;
				UnitSpec.uboone_e_dirt[0]= 47*4/5;
				UnitSpec.uboone_e_dirt[1]=47*1/5;
				UnitSpec.icarus_e_dirt[0]= 67*4/5;
				UnitSpec.icarus_e_dirt[1]=67*1/5;
				UnitSpec.sbnd_e_cosmo[0] = 9  ;
				UnitSpec.uboone_e_cosmo[0]= 11;
				UnitSpec.icarus_e_cosmo[0]= 10;

			


				std::vector<double > pred6 = UnitSpec.get_sixvector();
				std::vector<double > pred9 = UnitSpec.get_ninevector();
				std::vector<double > pred = UnitSpec.get_vector();

			

				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
				double chi2=0;
			
				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();


				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += (back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				std::cout<<unitModel.UUee<<" "<<unitModel.UUmm<<" "<<unitModel.UUem<<" "<<unitModel.UUme<<" "<<chi2<<std::endl;


		}}}


}














if(fraction_flag)
{
	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;

	if(dis_flag && !app_flag){
		usedetsys=false;
	}

	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
	bkgspec.sbnd_e_dirt[1]=44*1/5;
	bkgspec.uboone_e_dirt[0]= 47*4/5;
	bkgspec.uboone_e_dirt[1]=47*1/5;
	bkgspec.icarus_e_dirt[0]= 67*4/5;
	bkgspec.icarus_e_dirt[1]=67*1/5;
	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;
	

	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();

		int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */
		TMatrixT <double> M(matrix_size,matrix_size);
		TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
		TMatrixT <double> McI(matrix_size_c, matrix_size_c);


		int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int contMsize = (N_e_bins+N_m_bins)*N_dets;

		TMatrixT <double> Msys(bigMsize,bigMsize);
		sys_fill(Msys,usedetsys);

		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back[i]*back[j];
			}
		}




		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back);

		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}
		//Mtotal = Mstat;

		TMatrixT<double > Mctotal(contMsize,contMsize);
		contract_signal2(Mtotal,Mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		McI = Mctotal.Invert(&invdet);

		std::vector<std::vector<double >> vMcI = to_vector(McI);

	char filename[200];
	if(num_ster == 1){
		sprintf(filename,"GlobalFits/ntuples/nt_31_all_processed.root"); 
	} else if (num_ster == 2){
		sprintf(filename,"GlobalFits/ntuples/nt_32_all_processed.root"); 

	} else if(num_ster == 3){
		sprintf(filename,"GlobalFits/ntuples/nt_33_all_processed.root"); 

	}
	char outfilename[200];
	sprintf(outfilename,"ntuples/nt_3%d_all_processed_SBN.root",num_ster);

		TFile outputFile(outfilename,"RECREATE");
//		outputFile.cd();
		TNtuple ntuple("SBN1_99","SBN1_99","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56:mychi2");



	//std::cout<<filename<<std::endl;
	TFile *fm= new TFile(filename);
	TTree *chi2_90 =(TTree*)fm->Get("chi2_99_pr");
	 Float_t chi2, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45,phi46,phi56;
	 Float_t m4 = 0;
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
	 int nentries = chi2_90->GetEntries();
	 for (int i=0;i<nentries;i++) {
	        chi2_90->GetEntry(i);
	
	//	std::cout<<i<<" input_mn: "<<m4<<" "<<m5<<" "<<m6<<" input_ue "<<ue4<<" "<<ue5<<" "<<ue6<<" input_um4: "<<um4<<" "<<um5<<" "<<um6<<" input_chi: "<<chi2<<" "<<std::endl;

				double imn[3] = {(double)m4,(double)m5,(double)m6};
				double iue[3] = {ue4,ue5,ue6};
				double ium[3] = {um4, um5, um6};
				double iph[3] = {phi45,phi46, phi45};

				neutrinoModel appearanceModel(imn,iue,ium,iph);
				
				SBN_spectrum AppSpec(appearanceModel);
				//std::cout<<AppSpec.workingModel.mNu[0]<<" "<<AppSpec.workingModel.mNu[1]<<" "<<AppSpec.workingModel.mNu[2]<<std::endl;

				AppSpec.load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
				AppSpec.load_freq_3p3(SBND);
				AppSpec.load_freq_3p3(UBOONE);

					
				AppSpec.sbnd_e_dirt[1]=44*1/5;
				AppSpec.uboone_e_dirt[0]= 47*4/5;
				AppSpec.uboone_e_dirt[1]=47*1/5;
				AppSpec.icarus_e_dirt[0]= 67*4/5;
				AppSpec.icarus_e_dirt[1]=67*1/5;
				AppSpec.sbnd_e_cosmo[0] = 9  ;
				AppSpec.uboone_e_cosmo[0]= 11;
				AppSpec.icarus_e_cosmo[0]= 10;

			


				std::vector<double > pred6 = AppSpec.get_sixvector();
				std::vector<double > pred9 = AppSpec.get_ninevector();
				std::vector<double > pred = AppSpec.get_vector();

			

				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
				double mychi2=0;
			
				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						mychi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}


			std::cout<<i<<" mn: "<<AppSpec.workingModel.mNu[0]<<" "<<AppSpec.workingModel.mNu[1]<<" "<<AppSpec.workingModel.mNu[2]<<" electron: "<<ue4<<" "<<ue5<<" "<<ue6<<" muon: "<<um4<<" "<<um5<<" "<<um6<<" phi: "<<phi45<<" "<<phi46<<" "<<phi56<<" intpu_chi: "<<chi2<<" output_chi: "<<mychi2<<" "<<std::endl;
			ntuple.Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,mychi2);	
	
	
	
	
	
	
	
	
	
	
	
	
	 }
	fm->Close();
		outputFile.cd();
		ntuple.Write();
   		outputFile.Close();

}

if(fraction_flag&& false)
{
	std::cout<<"Starting fraction POT run"<<std::endl;
	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;

	if(dis_flag && !app_flag){
		usedetsys=false;
	}
	int ncovered = 0;
	int ncovered2 = 0;
	double FRAC_SIG_3 = 9.21+0.28;
	double FRAC_SIG_2 = 4.61 +0.28;
	double FRAC_SIG_1 = 1+0.28;

	char filename[200];
	if(num_ster == 1){
		sprintf(filename,"GlobalFits/ntuples/nt_31_all_processed.root"); 
	} else if (num_ster == 2){
		sprintf(filename,"GlobalFits/ntuples/nt_32_all_processed.root"); 

	} else if(num_ster == 3){
		sprintf(filename,"GlobalFits/ntuples/nt_33_all_processed.root"); 

	}
	char outfilename[200];
	sprintf(outfilename,"ntuples/nt_3%d_all_processed_SBN.root",num_ster);

		//TFile outputFile(outfilename,"RECREATE");
//		outputFile.cd();
		TNtuple ntuple("SBN1_99","SBN1_99","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56:mychi2");



	//std::cout<<filename<<std::endl;
	TFile *fm= new TFile(filename);
	TTree *chi2_90 =(TTree*)fm->Get("chi2_99_pr");
	 Float_t chi2, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45,phi46,phi56;
	 Float_t m4 = 0;
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
	 
	 int nentries = chi2_90->GetEntries();


	double pot = 0;
	for(double ipot = -5; ipot<0.25;ipot+=0.25){
	 ncovered = 0;
	 ncovered2=0;

	 pot = pow(10,ipot);
 


	 neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
	bkgspec.sbnd_e_dirt[1]=44*1/5;
	bkgspec.uboone_e_dirt[0]= 47*4/5;
	bkgspec.uboone_e_dirt[1]=47*1/5;
	bkgspec.icarus_e_dirt[0]= 67*4/5;
	bkgspec.icarus_e_dirt[1]=67*1/5;
	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;
	
				for(int i = 0; i < N_e_bins; i++){

					
					bkgspec.sbnd_e[i]= bkgspec.sbnd_e[i]*pot;
					bkgspec.sbnd_e_pho[i]= bkgspec.sbnd_e_pho[i]*pot;
					bkgspec.sbnd_e_dirt[i]= bkgspec.sbnd_e_dirt[i]*pot;
					bkgspec.sbnd_e_mu[i]= bkgspec.sbnd_e_mu[i]*pot;
					bkgspec.sbnd_f[i]= bkgspec.sbnd_f[i]*pot;
					bkgspec.sbnd_f_bar[i]= bkgspec.sbnd_f_bar[i]*pot;
	
					bkgspec.uboone_e[i]= bkgspec.uboone_e[i]*(0.5+pot);
					bkgspec.uboone_e_pho[i]= bkgspec.uboone_e_pho[i]*(0.5+pot);
					bkgspec.uboone_e_dirt[i]= bkgspec.uboone_e_dirt[i]*(0.5+pot);
					bkgspec.uboone_e_mu[i]= bkgspec.uboone_e_mu[i]*(0.5+pot);
					bkgspec.uboone_f[i]= bkgspec.uboone_f[i]*(0.5+pot);
					bkgspec.uboone_f_bar[i]= bkgspec.uboone_f_bar[i]*(0.5+pot);			
				
					bkgspec.icarus_e[i]= bkgspec.icarus_e[i]*pot;
					bkgspec.icarus_e_pho[i]= bkgspec.icarus_e_pho[i]*pot;
					bkgspec.icarus_e_dirt[i]= bkgspec.icarus_e_dirt[i]*pot;
					bkgspec.icarus_e_mu[i]= bkgspec.icarus_e_mu[i]*pot;
					bkgspec.icarus_f[i]= bkgspec.icarus_f[i]*pot;
					bkgspec.icarus_f_bar[i]= bkgspec.icarus_f_bar[i]*pot;
				
				}
			
				for(int i =0; i< N_m_bins; i++){

					bkgspec.sbnd_m[i]= bkgspec.sbnd_m[i]*pot;
					bkgspec.sbnd_m_pion[i]= bkgspec.sbnd_m_pion[i]*pot;

					bkgspec.uboone_m[i]= bkgspec.uboone_m[i]*(0.5+pot);
					bkgspec.uboone_m_pion[i]= bkgspec.uboone_m_pion[i]*(0.5+pot);

					bkgspec.icarus_m[i]= bkgspec.icarus_m[i]*pot;
					bkgspec.icarus_m_pion[i]= bkgspec.icarus_m_pion[i]*pot;

				}
	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();

		int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */
		TMatrixT <double> M(matrix_size,matrix_size);
		TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
		TMatrixT <double> McI(matrix_size_c, matrix_size_c);


		int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int contMsize = (N_e_bins+N_m_bins)*N_dets;

		TMatrixT <double> Msys(bigMsize,bigMsize);
		sys_fill(Msys,usedetsys);

		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back[i]*back[j];
			}
		}




		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back);

		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}
		//Mtotal = Mstat;

		TMatrixT<double > Mctotal(contMsize,contMsize);
		contract_signal2(Mtotal,Mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		McI = Mctotal.Invert(&invdet);

		std::vector<std::vector<double >> vMcI = to_vector(McI);




	ncovered = 0;
	ncovered2 = 0;
	 for (int i=0;i<nentries;i++) {
	        chi2_90->GetEntry(i);
			
	//	std::cout<<i<<" input_mn: "<<m4<<" "<<m5<<" "<<m6<<" input_ue "<<ue4<<" "<<ue5<<" "<<ue6<<" input_um4: "<<um4<<" "<<um5<<" "<<um6<<" input_chi: "<<chi2<<" "<<std::endl;

				double imn[3] = {(double)m4,(double)m5,(double)m6};
				double iue[3] = {ue4,ue5,ue6};
				double ium[3] = {um4, um5, um6};
				double iph[3] = {phi45,phi46, phi45};

				neutrinoModel appearanceModel(imn,iue,ium,iph);
				
				SBN_spectrum AppSpec(appearanceModel);
				//std::cout<<AppSpec.workingModel.mNu[0]<<" "<<AppSpec.workingModel.mNu[1]<<" "<<AppSpec.workingModel.mNu[2]<<std::endl;

				AppSpec.load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
				AppSpec.load_freq_3p3(SBND);
				AppSpec.load_freq_3p3(UBOONE);

					
				AppSpec.sbnd_e_dirt[1]=44*1/5;
				AppSpec.uboone_e_dirt[0]= 47*4/5;
				AppSpec.uboone_e_dirt[1]=47*1/5;
				AppSpec.icarus_e_dirt[0]= 67*4/5;
				AppSpec.icarus_e_dirt[1]=67*1/5;
				AppSpec.sbnd_e_cosmo[0] = 9  ;
				AppSpec.uboone_e_cosmo[0]= 11;
				AppSpec.icarus_e_cosmo[0]= 10;

		
				// Change PoT!!

				for(int i = 0; i < N_e_bins; i++){

					
					AppSpec.sbnd_e[i]= AppSpec.sbnd_e[i]*pot;
					AppSpec.sbnd_e_pho[i]= AppSpec.sbnd_e_pho[i]*pot;
					AppSpec.sbnd_e_dirt[i]= AppSpec.sbnd_e_dirt[i]*pot;
					AppSpec.sbnd_e_mu[i]= AppSpec.sbnd_e_mu[i]*pot;
					AppSpec.sbnd_f[i]= AppSpec.sbnd_f[i]*pot;
					AppSpec.sbnd_f_bar[i]= AppSpec.sbnd_f_bar[i]*pot;
	
					AppSpec.uboone_e[i]= AppSpec.uboone_e[i]*(0.5+pot);
					AppSpec.uboone_e_pho[i]= AppSpec.uboone_e_pho[i]*(0.5+pot);
					AppSpec.uboone_e_dirt[i]= AppSpec.uboone_e_dirt[i]*(0.5+pot);
					AppSpec.uboone_e_mu[i]= AppSpec.uboone_e_mu[i]*(0.5+pot);
					AppSpec.uboone_f[i]= AppSpec.uboone_f[i]*(0.5+pot);
					AppSpec.uboone_f_bar[i]= AppSpec.uboone_f_bar[i]*(0.5+pot);			
				
					AppSpec.icarus_e[i]= AppSpec.icarus_e[i]*pot;
					AppSpec.icarus_e_pho[i]= AppSpec.icarus_e_pho[i]*pot;
					AppSpec.icarus_e_dirt[i]= AppSpec.icarus_e_dirt[i]*pot;
					AppSpec.icarus_e_mu[i]= AppSpec.icarus_e_mu[i]*pot;
					AppSpec.icarus_f[i]= AppSpec.icarus_f[i]*pot;
					AppSpec.icarus_f_bar[i]= AppSpec.icarus_f_bar[i]*pot;
				
				}
			
				for(int i =0; i< N_m_bins; i++){

					AppSpec.sbnd_m[i]= AppSpec.sbnd_m[i]*pot;
					AppSpec.sbnd_m_pion[i]= AppSpec.sbnd_m_pion[i]*pot;

					AppSpec.uboone_m[i]= AppSpec.uboone_m[i]*(0.5+pot);
					AppSpec.uboone_m_pion[i]= AppSpec.uboone_m_pion[i]*(0.5+pot);

					AppSpec.icarus_m[i]= AppSpec.icarus_m[i]*pot;
					AppSpec.icarus_m_pion[i]= AppSpec.icarus_m_pion[i]*pot;

				}




				std::vector<double > pred6 = AppSpec.get_sixvector();
				std::vector<double > pred9 = AppSpec.get_ninevector();
				std::vector<double > pred = AppSpec.get_vector();

			

				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
				double mychi2=0;
			
				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						mychi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}

				if(mychi2 > FRAC_SIG_3){
					ncovered++;
				}
				if(mychi2 > FRAC_SIG_2){
					ncovered2++;
				}

			//std::cout<<i<<" mn: "<<AppSpec.workingModel.mNu[0]<<" "<<AppSpec.workingModel.mNu[1]<<" "<<AppSpec.workingModel.mNu[2]<<" electron: "<<ue4<<" "<<ue5<<" "<<ue6<<" muon: "<<um4<<" "<<um5<<" "<<um6<<" phi: "<<phi45<<" "<<phi46<<" "<<phi56<<" intpu_chi: "<<chi2<<" output_chi: "<<mychi2<<" "<<std::endl;
		//	ntuple.Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,mychi2);	
	
	 }// end of above n_entries loop

	 std::cout<<pot<<" "<<ncovered<<" "<<nentries<<" "<<(double)ncovered/((double)nentries)<<" "<<ncovered2<<" "<<" "<<(double)ncovered2/((double)nentries)<<std::endl;
	




	}// end of pot loop
	fm->Close();
	//	outputFile.cd();
	//	ntuple.Write();
   	//	outputFile.Close();

}





//Begin program flow control
if(fit_flag){

/*************************************************************
 *************************************************************
 *		Begin SBN Fitting Code
 ************************************************************
 ************************************************************/

	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);

	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);

	bkgspec.sbnd_e_dirt[0] = 44  ;
	bkgspec.uboone_e_dirt[0]= 47;
	bkgspec.icarus_e_dirt[0]= 67;

	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;


	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();



				neutrinoModel wrkModel(sqrt(in_dm),in_ue4,in_um4);
				wrkModel.dm41Sq = in_dm;

				SBN_spectrum wrkSpec(wrkModel);

				wrkSpec.load_freq_3p3(ICARUS);	
				wrkSpec.load_freq_3p3(UBOONE);	
				wrkSpec.load_freq_3p3(SBND);

				wrkSpec.sbnd_e_dirt[0] = 44  ;
				wrkSpec.uboone_e_dirt[0]= 47;
				wrkSpec.icarus_e_dirt[0]= 67;
				wrkSpec.sbnd_e_cosmo[0] = 9;
				wrkSpec.uboone_e_cosmo[0]= 11;
				wrkSpec.icarus_e_cosmo[0]= 10;

				std::vector<double > pred6 = wrkSpec.get_sixvector();
				std::vector<double > pred9 = wrkSpec.get_ninevector();

				
				int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
				int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

				/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
				 * */
				TMatrixT <double> M(matrix_size,matrix_size);
				TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
				TMatrixT <double> McI(matrix_size_c, matrix_size_c);


				std::vector<double > pred = wrkSpec.get_vector();

				int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
				int contMsize = (N_e_bins+N_m_bins)*N_dets;


				TMatrixT <double> Msys(bigMsize, bigMsize);
				sys_fill(Msys,true);

				for(int i =0; i<Msys.GetNcols(); i++)
				{
					for(int j =0; j<Msys.GetNrows(); j++)
					{
						Msys(i,j)=Msys(i,j)*back[i]*back[j];
					}
				}




				TMatrixT <double> Mstat(bigMsize, bigMsize);
				stats_fill(Mstat, back);

				TMatrixT <double > Mtotal(bigMsize, bigMsize);
				Mtotal =Msys+Mstat;

				TMatrixT<double > Mctotal(contMsize,contMsize);
				contract_signal2(Mtotal,Mctotal);

				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
				double invdet=0; // just to hold determinant
				double chi2=0;
					
				//bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
				std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				std::cout<<wrkModel.dm41Sq<<" "<<wrkModel.Ue[0]<<" "<<wrkModel.Um[0]<<" "<<chi2<<" "<<std::endl;
		
				for(int i =0; i<back6.size(); i++){
				std::cout<<back6[i]<<" "<<pred6[i]<<" "<<pred6[i]-back6[i]<<std::endl;

				}


}// end fit flag


if(bkg_flag){
	/*************************************************************
	 *************************************************************
	 *			BKG running area,  
	 *		I.E reproducing the histograms in SBN proposal
	 *	   Currently cant plot. not sure why, so output .roots
	 ************************************************************
	 ************************************************************/
	system("rm bkg_data/ICARUS_-inf.root");
	system("rm bkg_data/SBND_-inf.root");
	system("rm bkg_data/uBooNE_-inf.root");
		
	// Create a model to test, as we are reproducing bkg, create NULL model (m=0, u=0)
	neutrinoModel nullModel(0,0,0);

	//and create a SBN spectrum with the model we want to test
	SBN_spectrum bkgspectrum(nullModel);

	/* bkgspectrum can then do lots of useful things
	 * We create three detectors objects, currently overloaded to be construct via integer identifier
	 * but can pass in baseline, width, mass ..etc for a custom detector (see detector.h) 
  	 *	SBN_detector * SBND = new SBN_detector(0);
	 *	SBN_detector * UBOONE = new SBN_detector(1);
	 * 	SBN_detector * ICARUS = new SBN_detector(2);
	 *
	 *we can then get the events for any single detector or channel via
	 * 	bkgspectrum.fill_app(SBN_detector *);
	 * or
	 * 	bkgspectrum.fill_dis(SBN_detector *);
	 */

	//Alternatively we can run oscillate()
	//will then run the model over all detectors
	// creating the detecter objects internally for all uses

	//bkgspectrum.oscillate();
	bkgspectrum.oscillate_sample();


	system("cp bkg_data/ICARUS_-inf.root bkg_data/ICARUS_bkg.root"); 
	system("cp bkg_data/SBND_-inf.root bkg_data/SBND_bkg.root"); 
	system("cp bkg_data/uBooNE_-inf.root bkg_data/uBooNE_bkg.root");

	//and keep the end histograms in internal variable, sbnd_e sbnd_m ..etc..	
	//as well as writing them all to root files for plotting.. etc..
	// NO lotting should be done in this sbnfit code. doesnt seem right. 

/*	std::cout<<"****************** mu-spectra *******************************"<<std::endl;
	for(int i = 0; i < bkgspectrum.sbnd_m.size(); i++)
	{
		std::cout<<bkgspectrum.sbnd_m[i]<<" "<<bkgspectrum.uboone_m[i]<<" "<<bkgspectrum.icarus_m[i]<<std::endl;
	}

	std::cout<<"****************** Intrinsic-spectra *******************************"<<std::endl;
	for(int i = 0; i < bkgspectrum.sbnd_e.size(); i++)
	{
		std::cout<<bkgspectrum.sbnd_e[i]<<" "<<bkgspectrum.uboone_e[i]<<" "<<bkgspectrum.icarus_e[i]<<std::endl;
	}

	std::cout<<"****************** full_osc-spectra *******************************"<<std::endl;
	for(int i = 0; i < bkgspectrum.sbnd_f.size(); i++)
	{
		std::cout<<bkgspectrum.sbnd_f[i]<<" "<<bkgspectrum.uboone_f[i]<<" "<<bkgspectrum.icarus_f[i]<<std::endl;
	}
*/
	/* can also produce a "9-vector", i.e a vector that contains the 
	 * three nu_f, nu_e and nu_mu spectra in a row (for three detectors) for the purposes of 
	 * a chi square say.
	 *  in this order sbnd_osc_e, uboone_osc_e, icarus_osc_e, sbnd_intrin_e, uboone_intrin_e, icarus_intrin_e, sbnd_mu, uboone_mu, icarus_mu
	 * */

	std::vector<double > prediction9 = bkgspectrum.get_ninevector();

	/* can also produce a "6-vector", i.e a vector that contains the 
	 * same information, collaposed to just e and mu (so e is intrinsic + oscillated)
	 * */

	std::vector<double > prediction6 = bkgspectrum.get_sixvector();




}// end bkg_flag

if(sample_flag)
{
	
	for(double m = -2.0; m <=2.04; m=m+0.04){
		std::cout<<"Starting run for DM41: "<<m<<std::endl;

		neutrinoModel sModel(sqrt(pow(10,m)),1.0,1.0);
		sModel.dm41Sq = pow(10,m);

		SBN_spectrum samplespectrum(sModel);


		samplespectrum.oscillate_sample();

	}

	

}

if(cov_flag){






}

if(sens_flag)
{

 	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;

	if(dis_flag && !app_flag){
		usedetsys=false;
	}

	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
	bkgspec.sbnd_e_dirt[1]=44*1/5;
	bkgspec.uboone_e_dirt[0]= 47*4/5;
	bkgspec.uboone_e_dirt[1]=47*1/5;
	bkgspec.icarus_e_dirt[0]= 67*4/5;
	bkgspec.icarus_e_dirt[1]=67*1/5;
	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;
	

	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();

		int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */
		TMatrixT <double> M(matrix_size,matrix_size);
		TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
		TMatrixT <double> McI(matrix_size_c, matrix_size_c);


		int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int contMsize = (N_e_bins+N_m_bins)*N_dets;

		TMatrixT <double> Msys(bigMsize,bigMsize);
		sys_fill(Msys,usedetsys);

		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back[i]*back[j];
			}
		}




		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back);

		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}
		//Mtotal = Mstat;

		TMatrixT<double > Mctotal(contMsize,contMsize);
		contract_signal2(Mtotal,Mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		McI = Mctotal.Invert(&invdet);
		std::vector<std::vector<double >> vMcI = to_vector(McI);



	std::cout<<"Begining sensitivty scan"<<std::endl;

	/*************************************************************
	 *************************************************************
	 *		Sensitivity Analysis, currently for n=1 
	 *************************************************************
	 ************************************************************/
	if(sens_num == 1 && app_flag)
	{

		for(double m = -2; m <=2.04; m=m+0.04){
			for(int i = 0; i< 500; i++){

				double uei = rangen->Uniform(-0.30103,-5.3);
				neutrinoModel appearanceModel(sqrt(pow(10,m)), pow(10,uei),1.0);

				appearanceModel.dm41Sq = pow(10,m);

				SBN_spectrum AppSpec(appearanceModel);
	
				AppSpec.load_freq(ICARUS,0);//0 is silly app flag (get rid of this)
				AppSpec.load_freq(SBND,0);
				AppSpec.load_freq(UBOONE,0);

					
				AppSpec.sbnd_e_dirt[1]=44*1/5;
				AppSpec.uboone_e_dirt[0]= 47*4/5;
				AppSpec.uboone_e_dirt[1]=47*1/5;
				AppSpec.icarus_e_dirt[0]= 67*4/5;
				AppSpec.icarus_e_dirt[1]=67*1/5;
				AppSpec.sbnd_e_cosmo[0] = 9  ;
				AppSpec.uboone_e_cosmo[0]= 11;
				AppSpec.icarus_e_cosmo[0]= 10;

			


				std::vector<double > pred6 = AppSpec.get_sixvector();
				std::vector<double > pred9 = AppSpec.get_ninevector();

				
			
				std::vector<double > pred = AppSpec.get_vector();

			

				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
				double chi2=0;
			
				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				double sin22em = 4.0*pow(appearanceModel.Ue[0]*appearanceModel.Um[0],2.0);

				std::cout<<m<<" "<<appearanceModel.Ue[0]<<" "<<appearanceModel.Um[0]<<" "<<chi2<<" "<<sin22em<<std::endl;
			}//end random u run
		}//end mass run
	} //end 3p1 APPearance only sensitivity analysis

	if(sens_num == 1 && dis_flag && dis_which==1)
	{
		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(int i = 0; i< 333; i++){

				double umi = rangen->Uniform(-0.14,-2.0);
				neutrinoModel disappearanceModel(sqrt(pow(10,m)), 0.0,pow(10,umi));

				disappearanceModel.dm41Sq = pow(10,m);


				SBN_spectrum DisSpec(disappearanceModel);
				
				DisSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				DisSpec.load_freq_3p3(SBND);
				DisSpec.load_freq_3p3(UBOONE);

				DisSpec.sbnd_e_dirt[0] = 44*4/5  ;
				DisSpec.sbnd_e_dirt[1]=44*1/5;
				DisSpec.uboone_e_dirt[0]= 47*4/5;
				DisSpec.uboone_e_dirt[1]=47*1/5;
				DisSpec.icarus_e_dirt[0]= 67*4/5;
				DisSpec.icarus_e_dirt[1]=67*1/5;
				DisSpec.sbnd_e_cosmo[0] = 9  ;
				DisSpec.uboone_e_cosmo[0]= 11;
				DisSpec.icarus_e_cosmo[0]= 10;



			

				std::vector<double > pred = DisSpec.get_vector();
				std::vector<double > pred6= DisSpec.get_sixvector();

			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}

					//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				double chi2=0;
				
				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				double sin22mm = 4.0*(1-pow(disappearanceModel.Um[0],2.0))*pow(disappearanceModel.Um[0],2.0);

				std::cout<<m<<" "<<disappearanceModel.Ue[0]<<" "<<disappearanceModel.Um[0]<<" "<<chi2<<" "<<sin22mm<<std::endl;
			}//end random um4
		}//end m for loop
	} //end 3p1 sensitivity disapearance only analysis 

	if(sens_num == 1 && dis_flag && dis_which == 0)// This is the nu_e disapearance only channel, Interesting
	{
		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(int i = 0; i< 333; i++){

				double uei = rangen->Uniform(-0.14,-2.0);
				neutrinoModel disappearanceModel(sqrt(pow(10,m)), pow(10,uei),0);

				disappearanceModel.dm41Sq = pow(10,m);
				SBN_spectrum DisSpec(disappearanceModel);
				
				DisSpec.load_freq_3p3(ICARUS);
				DisSpec.load_freq_3p3(SBND);
				DisSpec.load_freq_3p3(UBOONE);

				DisSpec.sbnd_e_dirt[0] = 44*4/5  ;
				DisSpec.sbnd_e_dirt[1]=44*1/5;
				DisSpec.uboone_e_dirt[0]= 47*4/5;
				DisSpec.uboone_e_dirt[1]=47*1/5;
				DisSpec.icarus_e_dirt[0]= 67*4/5;
				DisSpec.icarus_e_dirt[1]=67*1/5;
				DisSpec.sbnd_e_cosmo[0] = 9  ;
				DisSpec.uboone_e_cosmo[0]= 11;
				DisSpec.icarus_e_cosmo[0]= 10;

				std::vector<double > pred = DisSpec.get_vector();
				std::vector<double > pred6= DisSpec.get_sixvector();

			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}

					//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				double chi2=0;
				
				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				double sin22ee = 4.0*(1-pow(disappearanceModel.Ue[0],2.0))*pow(disappearanceModel.Ue[0],2.0);

				std::cout<<m<<" "<<disappearanceModel.Ue[0]<<" "<<chi2<<" "<<sin22ee<<std::endl;
			}//end random ue4
		}//end m for loop
	} //end 3p1 sensitivity disapearance (ue4 dis only) only analysis 




	if(sens_num == 1 && both_flag)
	{
	std::cout<<"Begining N=1, dual appearance and dissapearance"<<std::endl;

		std::cout<<"Initialising output files"<<std::endl;
		TFile outputFile("hmm.root","RECREATE");
		outputFile.cd();
		std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
		TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double m = -2.00; m <=2.04; m=m+0.08){
			for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
				
				neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				bothModel.dm41Sq = pow(10,m);
			
				
			//	double imn[3] = {sqrt(pow(10,m)),sqrt(pow(10,1.24)),0.0};
			//	double iue[3] = {uei,0.069,0};
			//	double ium[3] = {umi, 0.16, 0.0};
			//	double iph[3] = {1.8*3.14159,0.0, 0.0};

			//	neutrinoModel bothModel(imn,iue,ium,iph);
				


				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}

	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/

				ntuple.Fill(bothModel.dm41Sq,bothModel.Ue[0],bothModel.Um[0],chi2);
				std::cout<<bothModel.dm41Sq<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<std::endl;
			}//end um4
			std::cout<<"#Finished m: "<<m<<" "<<umi<<std::endl;
			}//end ue4
		}//end m for loop

 
		
		outputFile.cd();
		std::cout<<"write ntuple"<<std::endl;
		ntuple.Write();
	 	std::cout<<"close file"<<std::endl;
   		outputFile.Close();
		std::cout<<"end all"<<std::endl;
	} //end 3p1 sensitivity both analysis


if(sens_num == 2&& false)
	{
	std::cout<<"Begining N=2, dual appearance and dissapearance"<<std::endl;



		std::cout<<"Initialising output files"<<std::endl;
		TFile outputFile("hmm.root","RECREATE");
		outputFile.cd();
		std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
		TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double m = -2.00; m <=2.04; m=m+0.08){
			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			
				for(double sinsq = -4; sinsq <=0; sinsq = sinsq + 0.5){
				double umiMin = -10;
				double ueiMin = -10;
				double chiMin = 100000;
				double sinsqMin = -10;

			  

				for(int n = 0; n< 200; n++){

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
				
				double umi = rangen->Uniform(0.0,-4.0);

				double uei = log10(sqrt( pow(10,sinsq)/(4.0*pow(pow(10,umi),2)) )  );

				while( pow(10,uei)>1 ){
					umi = rangen->Uniform(0.0,-4.0);
					uei = log10(sqrt( pow(10,sinsq)/(4.0*pow(pow(10,umi),2) )  ));
					//std::cout<<umi<<" "<<uei<<" "<<sinsq<<std::endl;
				}
				
				double imn[3] = {sqrt(pow(10,m)),sqrt(pow(10,1.24)),0.0};
				double iue[3] = {pow(10,uei),0.069,0};
				double ium[3] = {pow(10,umi), 0.16, 0.0};
				double iph[3] = {1.8*3.14159, 0.0, 0.0};

				neutrinoModel bothModel(imn,iue,ium,iph);
				


				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}

	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*McI(i,j)*(back6[j]-pred6[j]);
					}
				}


				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
				if(chi2<chiMin){
						chiMin=chi2;
						umiMin=umi;
						ueiMin=uei;
						sinsqMin=sinsq;	
				}


				}//end random
				ntuple.Fill(pow(10,m),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m)<<" "<<pow(10,ueiMin)<<" "<<pow(10,umiMin)<<" "<<pow(10,sinsq)<<" "<<chiMin<<std::endl;
			}//end sinsq loop

				std::cout<<"#Finished m: "<<m<<" "<<std::endl;
		}//end m for loop

 
		
		outputFile.cd();
		std::cout<<"write ntuple"<<std::endl;
		ntuple.Write();
	 	std::cout<<"close file"<<std::endl;
   		outputFile.Close();
		std::cout<<"end all"<<std::endl;
	} //end 3p1 sensitivity both analysis
	
if(sens_num == 2&& false)   // This is the m41 m51 fixed phi case for best fit
	{
	std::cout<<"Begining N=2, Dm41 V Dm51"<<std::endl;



		//std::cout<<"Initialising output files"<<std::endl;
	//	TFile outputFile("hmm.root","RECREATE");
	//	outputFile.cd();
	//	std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
	//	TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double m4 = -2.00; m4 <=2.04; m4=m4+0.04){
		for(double m5 = m4; m5 <=2.04; m5=m5+0.04){
			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
			
				double minPhi = 0;
				double minChi = 10000;	
				double iphi = 1.8;
				
				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),0.0};
				double iue[3] = {0.15,0.13,0};  // These are best fit!
				double ium[3] = {0.069,0.16, 0.0};
				//double iue[3] = {0.2,0.2,0};  // generic ones
				//double ium[3] = {0.2,0.2,0.0};


				double iph[3] = {iphi*3.14159, 0.0, 0.0};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =2;

				double round54 = round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;

				if(fabs(bothModel.dm54Sq) >= 100){ 
				//	std::cout<<"skipping this one 1:"<<std::endl;
						continue;
				}
				if(round54 > 2 ){ 
				//	std::cout<<"skipping this one 1: round54 "<<round54<<std::endl;
						continue;
				}

//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}
			
			

				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<chi2<<std::endl;
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
		}//end m5 loop
				//std::cout<<"#Finished m: "<<m4<<" "<<std::endl;
		}//end m for loop

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis


	
if(sens_num == 2 )   // This is the m41 m51 margined phi case for best fit
	{
	std::cout<<"Begining N=2, Dm41 V Dm51: amrgin"<<std::endl;



		//std::cout<<"Initialising output files"<<std::endl;
	//	TFile outputFile("hmm.root","RECREATE");
	//	outputFile.cd();
	//	std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
	//	TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double m4 = 0.0; m4 <=2.04; m4=m4+0.04){
		for(double m5 = 0.0; m5 <=2.04; m5=m5+0.04){
			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
			
				double minPhi = 0;
				double minChi = 10000;	
				for(double iphi = 0; iphi< 2; iphi+=0.05){
				//double iphi = 1.8;
				
				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),0.0};
				double iue[3] = {0.15,0.13,0};  // These are best fit!
				double ium[3] = {0.069,0.16, 0.0};
				//double iue[3] = {0.2,0.2,0};  // generic ones
				//double ium[3] = {0.2,0.2,0.0};


				double iph[3] = {iphi*3.14159, 0.0, 0.0};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =2;

				double round54 = round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;

				if(fabs(bothModel.dm54Sq) >= 100){ 
				//	std::cout<<"skipping this one 1:"<<std::endl;
						continue;
				}
				if(round54 > 2 ){ 
				//	std::cout<<"skipping this one 1: round54 "<<round54<<std::endl;
						continue;
				}

//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}
			
				if(chi2<=minChi){
					minPhi=iphi*3.14159;
					minChi=chi2;
				}

					} // end phi loop


				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<minChi<<" "<<minPhi<<std::endl;
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
		}//end m5 loop
				//std::cout<<"#Finished m: "<<m4<<" "<<std::endl;
		}//end m for loop

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis

if(sens_num == 2&& false)   // This is the m51 and phi plot for averaged else
	{
	std::cout<<"Begining N=2, Dm51 V sinsq5"<<std::endl;



		//std::cout<<"Initialising output files"<<std::endl;
	//	TFile outputFile("hmm.root","RECREATE");
	//	outputFile.cd();
	//	std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
	//	TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double sinsq = -3.6; sinsq <=-2; sinsq+=0.1){
		for(double m5 = -2.0; m5 <=2.04; m5=m5+0.04){
			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			

				double m4 = -1;

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
				double minPhi = 0;	
				double minUm5 = 0;
				double minUe5 = 0;
				double minChi = 10000;	
			
				for(int n =0; n< 50; n++){	
		
				double iphi= rangen->Uniform(0,2);	
				double um5 = rangen->Uniform(0.0,-4.0);

				double ue5 = log10(sqrt( pow(10,sinsq)/(4.0*pow(pow(10,um5),2)) )  );

				while( pow(10,ue5) > 1 ){
					um5 = rangen->Uniform(0.0,-4.0);
					ue5 = log10(sqrt( pow(10,sinsq)/(4.0*pow(pow(10,um5),2) )  ));
					//std::cout<<umi<<" "<<uei<<" "<<sinsq<<std::endl;
				}


				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),0.0};
				//double iue[3] = {0.15,0.13,0};  // These are best fit!
				//double ium[3] = {0.069,0.16, 0.0};
				double iue[3] = {0.1,pow(10,ue5),0};//ue5,0};  // generic ones
				double ium[3] = {0.1,pow(10,um5),0};//um5,0.0};


				double iph[3] = {iphi*3.14159, 0.0, 0.0};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =2;

				double round54 = round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;

				if(fabs(bothModel.dm54Sq) >= 100){ 
				//	std::cout<<"skipping this one 1:"<<std::endl;
						continue;
				}
				if(round54 > 2 ){ 
				//	std::cout<<"skipping this one 1: round54 "<<round54<<std::endl;
						continue;
				}

//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);


					}
				}
			
				if(chi2<=minChi){
					minPhi=iphi*3.14159;
					minUe5= ue5;
					minUm5 = um5;
					minChi=chi2;

				}

				} // end N loop


				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<minChi<<" "<<minPhi<<" "<<minUe5<<" "<<minUm5<<" "<<sinsq<<std::endl;
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
		}//end m5 loop
				//std::cout<<"#Finished m: "<<m4<<" "<<std::endl;
		}//end sinsq

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis


if(sens_num == 2&&false)   // This is the m51 and phi plot for averaged else
	{
	std::cout<<"Begining N=2, Dm51 V phi"<<std::endl;



		//std::cout<<"Initialising output files"<<std::endl;
	//	TFile outputFile("hmm.root","RECREATE");
	//	outputFile.cd();
	//	std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
	//	TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double iphi =0; iphi <=2.025; iphi+=0.025){
		for(double m5 = -2.0; m5 <=2.04; m5=m5+0.04){
			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			

				double m4 = -1;

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
				double minPhi = 0;	
				double minUm5 = 0;
				double minUe5 = 0;
				double minChi = 10000;	
			
		
		

				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),0.0};
				//double iue[3] = {0.15,0.13,0};  // These are best fit!
				//double ium[3] = {0.069,0.16, 0.0};
				double iue[3] = {0.13,0.13,0};//ue5,0};  // generic ones
				double ium[3] = {0.05,0.05,0};//um5,0.0};


				double iph[3] = {iphi*3.14159, 0.0, 0.0};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =2;

				double round54 = round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;

				if(fabs(bothModel.dm54Sq) >= 100){ 
				//	std::cout<<"skipping this one 1:"<<std::endl;
						continue;
				}
				if(round54 > 2 ){ 
				//	std::cout<<"skipping this one 1: round54 "<<round54<<std::endl;
						continue;
				}

//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}


	
		
				double invdet=0; // just to hold determinant
				double chi2=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				double mod = 1.0;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);


					}
				}
			
			

				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<chi2<<" "<<iphi*3.14159<<std::endl;
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
		}//end m5 loop
				//std::cout<<"#Finished m: "<<m4<<" "<<std::endl;
		}//end sinsq

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis




int nn = 0;
if(sens_num == 3)
	{
	std::cout<<"Begining N=3, Dm41 V Dm51"<<std::endl;



		//std::cout<<"Initialising output files"<<std::endl;
	//	TFile outputFile("hmm.root","RECREATE");
	//	outputFile.cd();
	//	std::cout<<"Initialising ntuple ~ TNuple"<<std::endl;
	//	TNtuple ntuple("3p1chiNtuple","3p1chiNtuple","logDm4:logUe4:logUm4:Chi2");

                           /*     Double_t nUE4;
                                  Double_t nUM4;
                                  Double_t nDM4;
    				  Double_t nCHI;	

                                  ntuple->Branch("logUe4",&nUE4);
                                  ntuple->Branch("logUm4",&nUM4);
                                  ntuple->Branch("logDm4",&nDM4);
                                  ntuple->Branch("Chi",&nCHI);*/
     

		for(double m4 = -2.00; m4 <=2.04; m4+=0.04){
		for(double m5 = m4; m5 <=2.04; m5=m5+0.04){
	//	for(double m6 = m4; m6 <=2.04; m6+=0.04){

			nn++;
			//if(nn%2==0){continue;}

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
				
				double m6 = round(log10(22.0)/0.04)*0.04;
			//	double m5 = round(log10(17.0)/0.04)*0.04;
			

	//			for(double iphi = 0; iphi<2; iphi+=0.2){

				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),sqrt(pow(10,m6))};
				double iue[3] = {0.11,0.11,0.11};
				double ium[3] = {0.12,0.17, 0.14};
				double iph[3] = {1.6*3.14159, 0.28 *3.14159, 1.4*3.14159};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =3;

				double round54 = sgn(bothModel.dm54Sq)*round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;
				double round64 = sgn(bothModel.dm64Sq)*round(log10(fabs(bothModel.dm64Sq))/0.04)*0.04;
				double round65 = sgn(bothModel.dm65Sq)*round(log10(fabs(bothModel.dm65Sq))/0.04)*0.04;

				if(fabs(bothModel.dm54Sq) >= 100 || round54 > 2 ){
						continue;
				}
				if(fabs(bothModel.dm64Sq) >= 100 || round64 > 2 ){
						continue;
				}
				if(fabs(bothModel.dm65Sq) >= 100 || round65 > 2 ){
						continue;
				}
//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothSpec.sbnd_e_dirt[0] = 44*4/5  ;
				BothSpec.sbnd_e_dirt[1]=44*1/5;
				BothSpec.uboone_e_dirt[0]= 47*4/5;
				BothSpec.uboone_e_dirt[1]=47*1/5;
				BothSpec.icarus_e_dirt[0]= 67*4/5;
				BothSpec.icarus_e_dirt[1]=67*1/5;
				BothSpec.sbnd_e_cosmo[0] = 9 ;
				BothSpec.uboone_e_cosmo[0]= 11;
				BothSpec.icarus_e_cosmo[0]= 10;


				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > pred9 = BothSpec.get_ninevector();

				
				std::vector<double > pred = BothSpec.get_vector();
			
				if(pred6.size()!=Mctotal.GetNcols()){std::cout<<"ERROR"<<std::endl;}

	
		
				double invdet=0; // just to hold determinant
				double chi2b=0;
		
				//	bit o inverting, root tmatrix seems perfectly fast	
				McI = Mctotal.Invert(&invdet);

				//check for previous known bug!
				if(false && matrix_size_c != pred6.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred6.size()<<" back: "<<back6.size()<<std::endl;	
				}

				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();

				if(whatsize != back6.size() || whatsize != pred6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
				}

				//	std::cout<<"test: "<<McI[4][5]<<" "<<vMcI[4][5]<<std::endl;




				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2b += (back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}
				double sum_of_elems=0;
				for (double n : pred)
	    			sum_of_elems += n;

				double sum_of_back=0;
				for (double n : back)
	    			sum_of_back += n;
				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<pow(10,m6)<<" "<<chi2b<<std::endl;
				//std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<pow(10,m6)<<" "<<chi2b<<" "<<sum_of_elems<<" "<<sum_of_back<<" "<<sum_of_back-sum_of_elems<<" "<<m4<<" "<<m5<<std::endl;
//		}//phi loop
			
	
		}//end m5 loop
				//std::cout<<"#Finished m: "<<m4<<" "<<std::endl;
		}//end m for loop

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis
}//end sens_flag

	








if(test_flag){
	/*************************************************************
	 *************************************************************
	 *		test flag please ignore, 
	 *	no really, just a jumble place for running tests
	 ************************************************************
	 ************************************************************/


	std::cout<<"Running test Mode"<<std::endl;

 	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);


	bkgspec.sbnd_e_dirt[0] = 44  ;
	bkgspec.uboone_e_dirt[0]= 47;
	bkgspec.icarus_e_dirt[0]= 67;

	bkgspec.sbnd_e_cosmo[0] = 9  ;
	bkgspec.uboone_e_cosmo[0]= 11;
	bkgspec.icarus_e_cosmo[0]= 10;


	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();
	TRandom *rangen    = new TRandom();






if(true){

	neutrinoModel testModel(sqrt(pow(10,1.72)),100*pow(10,-4.95174),1);
	testModel.dm41Sq = pow(10,1.72);
	SBN_spectrum myspec(testModel);


	myspec.load_freq(ICARUS,0);
	myspec.load_freq(SBND,0);
	myspec.load_freq(UBOONE,0);


	std::vector<double > pred9 = myspec.get_vector();
	
	neutrinoModel testModelb(sqrt(0.11),25*pow(10,-2.64916),1);
	testModelb.dm41Sq = 0.1;
	SBN_spectrum myspecb(testModelb);


	//myspecb.load_freq(ICARUS,1);
	//myspecb.load_freq(SBND,1);
	//myspecb.load_freq(UBOONE,1);

	myspecb.load_freq_3p3(ICARUS);


	std::vector<double > pred9b = myspecb.get_vector();



	for(int i=0; i<pred9b.size(); i++)
	{

	std::cout<<back9[i]<<" "<<pred9[i]<<" "<<pred9b[i]<<std::endl;

	}

}


/*

	double Sapp=4*pow(testModel.Ue[0]*testModel.Um[0],2.0);
	double Sdis=1.0-4*pow(testModel.Um[0],2)*(1- pow(testModel.Um[0],2));
	std::cout<<"Ptest2 "<<Pmue(600,1.0,1.0,Sdis)<<std::endl;


	myspec.oscillate();
	myspec.vec_print();



*/



/*	SBN_spectrum nullsp;
	nullsp.vec_print();	
	std::cout<<"Begin Test_Flag"<<std::endl;
	

	std::cout<<"End SBN_spectrum Initilise"<<std::endl;
	neutrinoModel testModel;
	testModel.Ue[0]=0.0;
	testModel.Um[0]=0.0;
	testModel.mNu[0]=2.0;

	SBN_spectrum myspec(testModel);
	std::vector<double > pred = myspec.add_SBN_spectrum(nullsp);

		//myspec.vec_print();


	myspec.oscillate();

	std::vector<double > predtest = myspec.add_SBN_spectrum(nullsp);


	for(int i = 0; i<pred.size(); i++){
		std::cout<<predtest[i]<<std::endl;
	}
*/

} //end test flag





return 0;
}// end main





