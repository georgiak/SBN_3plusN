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
#include "TRandom.h"
#include "TError.h"

#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

#include "params.h"
#include "prob.h"
#include "model.h"
#include "correlation.h"
#include "wrkInstance.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


/*************************************************************
 *************************************************************
 *		Define a working instance class as I'm sick
 *		of re-nitialising and doing chi^2 every time
 ************************************************************
 ************************************************************/



/*************************************************************
 *************************************************************
 *		BEGIN Main::sbnfit.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


	/*	int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti;
		int contMsize = (N_e_bins+N_m_bins)*N_dets*N_anti;
		
		TMatrixT <double> Msys(bigMsize,bigMsize);
		TMatrixT <double> MsysC(contMsize,contMsize);

		for(int i =0; i< Msys.GetNcols(); i++){
		for(int j =0; j< Msys.GetNcols(); j++){
		Msys(i,j) = 1;
		}}

		contract_signal2_anti(Msys,MsysC);
		
	for(int i =0; i< MsysC.GetNcols(); i++){
		for(int j =0; j< MsysC.GetNcols(); j++){
			std::cout<<i<<" "<<j<<" "<<MsysC(i,j)<<std::endl;
		}}
exit(EXIT_FAILURE);
*/









//gSystem->Load("libTree");


// Just a few flags to control program flow.
bool fit_flag = false;
bool verbose_flag = false;
bool test_flag=false;
bool bkg_flag= false;
bool sample_flag = false;
bool cov_flag = false;
bool sens_flag=false;
bool comb_flag = false;

bool both_flag = true;
bool dis_flag = false;
bool app_flag = false;
bool wierd_flag = false;
int which_channel = BOTH_ONLY;

bool unit_flag = false;
bool fraction_flag = false;
bool anti_flag = false;
int anti_mode = 0;
bool inject_flag = false;


double bfmn[3] = {0.398107,1.0,0};
double bfue[3] = {0.13,0.14,0};
double bfum[3] = {0.15,0.13,0};
double bfphi[3] = {0.0,0.0,0.0};


neutrinoModel inputModel(bfmn,bfue,bfum,bfphi);

bool margin = false;


int plotmode = 1;

bool pot_flag = false;

int mode_flag = 0;

int parallel_split=0;

int  sens_num=1;
double dm41 = -1.0;

double in_dm = 0;
double in_ue4 = 0;
double in_um4=0;


double pot_num =0;
double pot_num_bar =0;

double inPhi45 = 0;

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
	{"mode",		required_argument,	0, 'M'},
	{"test", 		no_argument, 		0, 'T'},
	{"mass",		required_argument,	0, 'm'},
	{"anti",		no_argument,		0, 'A'},
	{"ue4", 		required_argument,	0, 'e'},
	{"help",		no_argument, 		0, 'h'},
	{"um4"	,		required_argument,	0, 'u'},
	{"inject",	required_argument,		0, 'I'},
	{"verbose",		no_argument, 		0, 'v'},
	{"sensitivity",		required_argument,	0, 'S'},
	{"stat-only",		no_argument,		0, 'l'},
	{"sample",		no_argument,		0, 's'},
	{"cov",			no_argument, 		0, 'c'},
	{"dis",			no_argument,	 	0, 'd'},
	{"wierd",		no_argument,		0, 'w'},
	{"signal",		required_argument,	0, 'g'},
	{"both",		no_argument,		0, 'b'},
	{"unitary",		no_argument,		0, 'n'},
	{"num",			required_argument,	0, 'N'},
	{"fraction",		required_argument,	0, 'f'},
	{"pot",			required_argument,	0, 'p'},
	{"plotmode",		required_argument,	0, 'k'},
	{"margin",		no_argument,		0, 'r'},
	{"app",			no_argument,		0, 'a'},
	{"phi",			required_argument,	0, 'P'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "I:daP:lf:nuM:e:m:svp:hS:crFTABN:bwk:g:", longopts, &index);

	switch(iarg)
	{
		case 'r':
			margin=true;
			break;
		case 'g':
			{
			
			double tm[3] ={0,0,0};
			double te[3] ={0,0,0};
			double tu[3] ={0,0,0};
			double tp[3] ={0,0,0};

			std::string s = optarg;
			std::string delimiter = ":";
			int cnt=0;
			size_t pos = 0;
			std::string token;
			while ((pos = s.find(delimiter)) != std::string::npos) {
			    token = s.substr(0, pos);
			    //std::cout << token << std::endl;
			    s.erase(0, pos + delimiter.length());

			    if(cnt<3){
				double mm=pow(10,round(log10(fabs( atof(token.c_str())))/0.04)*0.04);	
				tm[cnt] =mm; 
			    }
 			    if(cnt>=3&&cnt<6){
				te[cnt-3] = atof(token.c_str());	
			    } 
 			    if(cnt>=6&&cnt<9){
				tu[cnt-6] = atof(token.c_str());	
			    } 
			    if(cnt>=9&&cnt<12){
				tp[cnt-9] = atof(token.c_str());	
			    } 
			   
			    cnt++;
			}
			//std::cout << s << std::endl;//last one remains in s
			tp[2]=atof(s.c_str());
					
			inputModel=neutrinoModel(tm,te,tu,tp);

			}
			break;
		case 'P':
			inPhi45 = strtof(optarg,NULL);
			break;
		case 'I':
			{			
			inject_flag = true;
			std::string s = optarg;
			std::string delimiter = ":";
			//std::cout<<s.find(delimiter)<<std::endl;
			std::string spot = s.substr(0, s.find(delimiter));	
			std::string spotbar = s.substr(s.find(delimiter)+1);

			pot_num= atof(spot.c_str());
			pot_num_bar = atof(spotbar.c_str());
			}
		

			break;
		case 'F':
			fit_flag = true;
			//mS = strtof(optarg,NULL);
			break;
		case 'B':
			bkg_flag = true;
			break;
		case 'M':
			mode_flag = strtof(optarg,NULL);
			break;
		case 'A':
			anti_flag = true;
			anti_mode = 1;
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
		case 'p':
			{
			pot_flag = true;
			std::string s = optarg;
			std::string delimiter = ":";
			//std::cout<<s.find(delimiter)<<std::endl;
			std::string spot = s.substr(0, s.find(delimiter));	
			std::string spotbar = s.substr(s.find(delimiter)+1);

			pot_num= atof(spot.c_str());
			pot_num_bar = atof(spotbar.c_str());
			}
			break;
		case 'k':
			plotmode = strtof(optarg,NULL);
			break;
		case 'N':
			num_ster = strtof(optarg,NULL);
			break;
		case 's':
			sample_flag = true;
			break;
		case 'l':
			stat_only = true;
			break;
		case 'd':
			dis_flag = true;
			//	dis_which =strtof(optarg,NULL); //obsolete!
			which_channel = DIS_ONLY;
			break;
		case 'a':
			app_flag = true;
			which_channel = APP_ONLY;
			break;
		case 'w':
			wierd_flag = true;
			which_channel = WIERD_ONLY;
			break;
		case 'b':
			both_flag = true;
			which_channel = BOTH_ONLY;
			break;
		case 'm':
			in_dm  = strtof(optarg,NULL);
			break;
		case 'f':
			fraction_flag = true;
			parallel_split  = strtof(optarg,NULL);
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
			std::cout<<"\t\t\t --sensitivity 1 --app --dis 1 --mode 1,  makes 3p1_both_em.dat"<<std::endl;
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
			std::cout<<"\t-A\t--anti\t\t Anti neutrino running mode"<<std::endl;
			return 0;
	}

}



if(verbose_flag)
{
	std::cout<<"#v:***************** Initial Parameters ******************"<<std::endl;
	std::cout<<"#v: N_m_bins "<<N_m_bins<<std::endl;
	std::cout<<"#v: N_e_bins "<<N_e_bins<<std::endl;
	std::cout<<"#v: N_dets "<<N_dets<<std::endl;
	std::cout<<"#v: N_e_spectra "<<N_e_spectra<<std::endl;
	std::cout<<"#v: N_m_spectra "<<N_m_spectra<<std::endl;
	std::cout<<"#v: N_anti "<<N_anti<<std::endl;
	std::cout<<"#v: num sterile "<<num_ster<<std::endl;
	if(!anti_flag){
	std::cout<<"#v: NEUTRINO mode only "<<num_ster<<std::endl;
	}else{
	std::cout<<"#v: NU+NUBAR modes "<<num_ster<<std::endl;
	}	
	switch(which_channel){
			case APP_ONLY:
				std::cout<<"#v: Channel app";
				break;
			case DIS_ONLY:
				std::cout<<"#v: Channel dis";
				break;
			case BOTH_ONLY:
				std::cout<<"#v: Channel both";
				break;
			case WIERD_ONLY:
				std::cout<<"#v: Channel wierd";
				break;
	}
		
	std::cout<<"#*******************************************************"<<std::endl;
}




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










if(fraction_flag && !false ) // This is just an obsolete old one for reading it in and doing noting
{
	std::cout<<"filename"<<std::endl;
	char filename[200];
	if(num_ster == 1){
		sprintf(filename,"GlobalFits/ntuples/nt_31_all_processed.root");
//	        sprintf(filename,"GlobalFits/ntuples/nt_31_all.root");	
	} else if (num_ster == 2){

		sprintf(filename,"GlobalFits/ntuples/nt_32_all_processed.root"); 
	//	sprintf(filename,"GlobalFits/ntuples/nt_32_all.root"); 
	} else if(num_ster == 3){
		sprintf(filename,"GlobalFits/ntuples/nt_33_all_processed.root"); 

	}


	std::cout<<"read in`"<<std::endl;
	std::cout<<filename<<std::endl;
	TFile *fm= new TFile(filename);
//	TTree *chi2_90 =(TTree*)fm->Get("chi2_99");
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
	double mymax = 0;
	double mu[6] = {0,0,0,0,0,0};
 for (int i=0;i<nentries;i++) {
	        chi2_90->GetEntry(i);
	double tt = ue4*ue5*um4*um5;
		if(tt > mymax){
			mymax = tt;
			mu[0]=ue4;
			mu[1]=ue5;
			mu[2]=um4;
			mu[3]=um5;
			mu[4]=m4;
			mu[5]=m5;
			 std::cout<<"temp"<<mymax<<" "<<mu[0]<<" "<<mu[1]<<" "<<mu[2]<<" "<<mu[3]<<" "<<mu[4]<<" "<<mu[5]<<std::endl;
			 std::cout<<"dm41^2 : "<<mu[4]*mu[4]<<" dm51^2: "<<mu[5]*mu[5]<<" dm54^2: "<<mu[5]*mu[5]-mu[4]*mu[4]<<std::endl;
		

		}

		//std::cout<<m4<<" "<<sins2<<std::endl;
		//		//std::cout<<m4<<" "<<m5<<" "<<ue4<<" "<<ue5<<" "<<um4<<" "<<um5<<" "<<phi45<<" "<<chi2<<std::endl;
		//
					}
			 std::cout<<"biggest"<<mymax<<" "<<mu[0]<<" "<<mu[1]<<" "<<mu[2]<<" "<<mu[3]<<" "<<mu[4]<<" "<<mu[5]<<std::endl;
			 std::cout<<"dm41^2 : "<<mu[4]*mu[4]<<" dm51^2: "<<mu[5]*mu[5]<<std::endl;
			 
	return 0;	
	}


if(fraction_flag && !true) //this i smain!!
{

	double norm_pot = pow(10,pot_num);
	double norm_pot_bar = pow(10,pot_num_bar);// 1.0;//0.0000001;
	 
        wrkInstance fractionInstance(which_channel, anti_mode, norm_pot, norm_pot_bar);

//	exit(EXIT_FAILURE);

	char filename[200];
	if(num_ster == 1){
		sprintf(filename,"GlobalFits/ntuples/nt_31_all_processed.root");
	//        sprintf(filename,"GlobalFits/ntuples/nt_31_brute.root");	
	} else if (num_ster == 2){
		sprintf(filename,"GlobalFits/ntuples/nt_32_all_processed.root"); 

	} else if(num_ster == 3){
		sprintf(filename,"GlobalFits/ntuples/nt_33_all_processed.root"); 

	}

/*	char outfilename[200];
	sprintf(outfilename,"ntuples/nt_3%d_all_processed_SBN2.root",num_ster);

		TFile outputFile(outfilename,"RECREATE");
//		outputFile.cd();
		TNtuple ntuple("SBN1_99","SBN1_99","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56:mychi2");
*/

	TFile *fm= new TFile(filename);
	TTree *chi2_99 =(TTree*)fm->Get("chi2_99_pr");
	 Float_t chi2, ue4, um4, m5, ue5, um5, m6, ue6, um6, phi45,phi46,phi56;
	 Float_t m4 = 0;
         chi2_99->SetBranchAddress("chi2",&chi2);
         chi2_99->SetBranchAddress("m4",&m4);
         chi2_99->SetBranchAddress("ue4",&ue4);
         chi2_99->SetBranchAddress("um4",&um4);
         chi2_99->SetBranchAddress("m5",&m5);
         chi2_99->SetBranchAddress("ue5",&ue5);
         chi2_99->SetBranchAddress("um5",&um5);
         chi2_99->SetBranchAddress("m6",&m6);
         chi2_99->SetBranchAddress("ue6",&ue6);
         chi2_99->SetBranchAddress("um6",&um6);
         chi2_99->SetBranchAddress("phi45",&phi45);
         chi2_99->SetBranchAddress("phi46",&phi46);
         chi2_99->SetBranchAddress("phi56",&phi56);
	 int nentries = chi2_99->GetEntries();

	 double imin=0;
	 double imax=nentries;


	//Rudamentary parallalisation scheme
	 int pfrac=nentries/8;
	 if(parallel_split!=0){
		switch(parallel_split){
			case 1:
				imin=0;
				imax=pfrac;
				break;
			case 2:
				imin=pfrac;
				imax=2*pfrac;
				break;
			case 3:
				imin=2*pfrac;
				imax=3*pfrac;
				break;
			case 4:
				imin=3*pfrac;
				imax=4*pfrac;
				break;
			case 5:
				imin=4*pfrac;
				imax=5*pfrac;
				break;
			case 6:
				imin=5*pfrac;
				imax=6*pfrac;
				break;
			case 7:
				imin=6*pfrac;
				imax=7*pfrac;
				break;
			case 8:
				imin=7*pfrac;
				imax=nentries;
				break;
			}
	}


	 for (int i=imin;i<imax;i++) {
	        chi2_99->GetEntry(i);
	
				double imn[3] = {(double)m4,(double)m5,(double)m6};
				double iue[3] = {ue4,ue5,ue6};
				double ium[3] = {um4, um5, um6};
				double iph[3] = {phi45,phi46,phi56};

				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile=num_ster;

				fractionInstance.calc_chi(signalModel, i, norm_pot, norm_pot_bar);
			
				//exit(EXIT_FAILURE);	
	 }


	fm->Close();
	/*	outputFile.cd();
		ntuple.Write();
   		outputFile.Close();
	*/

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
		//sprintf(filename,"GlobalFits/ntuples/nt_31_all_processed.root"); 
 	        sprintf(filename,"GlobalFits/ntuples/nt_31_brute.root");	
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
	//TTree *chi2_90 =(TTree*)fm->Get("chi2_99_pr");
	TTree *chi2_90 =(TTree*)fm->Get("chi2_99");
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
			
	bkgspec.scale_by_pot(pot);


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



		
				// Change PoT!!
				AppSpec.scale_by_pot(pot);


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

if(pot_flag){

	double ipot = pow(10,pot_num);
	double ipotbar = pow(10,pot_num_bar);


	wrkInstance potInstance(which_channel, anti_mode, ipot, ipotbar);

	int vector_modifier = 1;
	
	if(anti_flag){
		vector_modifier =2;
	}


	//now load file with all other vectors.
	std::string filename = "3p1_both";

	if(num_ster == 1){
		switch(which_channel){
			case APP_ONLY:
				filename = "3p1_app";
				break;
			case DIS_ONLY:
				filename = "3p1_dis";
				break;
			case BOTH_ONLY:
				filename = "3p1_both";
				break;
			case WIERD_ONLY:
				filename = "3p1_wierd";
				break;
		}

	} else if (num_ster == 2){
			switch(which_channel){
			case APP_ONLY:
				filename = "3p2_app";
				break;
			case DIS_ONLY:
				filename = "3p2_dis";
				break;
			case BOTH_ONLY:
				filename = "3p2_both";
				break;
			case WIERD_ONLY:
				filename = "3p2_wierd";
				break;
		}
	} else if(num_ster == 3){
			switch(which_channel){
			case APP_ONLY:
				filename = "3p3_app";
				break;
			case DIS_ONLY:
				filename = "3p3_dis";
				break;
			case BOTH_ONLY:
				filename = "3p3_both";
				break;
			case WIERD_ONLY:
				filename = "3p3_wierd";
				break;
		}

	}

	std::string pre =  "fractiondata/";
	std::string post = ".dat";

	if(anti_flag){
		pre = "fractiondata/NUBAR_MODE/";	
		post = ".bar.dat";
	}
	pre.append(filename);
	pre.append(post);
	filename=pre;


if(filename != "none"){

	int k = 0;
	std::string num;
	std::string m4;
	std::string m5;
	std::string m6;

	std::string ue4;
	std::string ue5;
	std::string ue6;

	std::string um4;
	std::string um5;
	std::string um6;

	std::string phi45;
	std::string phi46;
	std::string phi56;

	std::string whicho;

	std::string pot;
	std::string oldchi;
	std::string newchi;

	std::string tempVec;	

	//no longer the pred6in, sometimes pred12 all or whatever bull
	std::vector<double > vectorin;

	double mn[3];
	double ue[3];
	double um[3];
	double phi[3];

	double potin;
	double oldchiin;
	double newchiin;

	int whichflagin = -1;

	int count = 0;

	std::ifstream myfile (filename);
	if (!myfile.is_open())
	{
		std::cout<<"#ERROR: Passed POT file does not exist: "<<filename<<std::endl;
		std::cout<<" which_channel "<<which_channel<<" numster "<<num_ster<<std::endl;
		exit(EXIT_FAILURE);
	}

		while(!myfile.eof()){
		
			count++;	
			myfile >> num;
			//std::cout<<num<<std::endl;	
		
			myfile >> m4;
			myfile >> m5;
			myfile >> m6;
			
			mn[0]= atof(m4.c_str());
			mn[1]= atof(m5.c_str());
			mn[2]= atof(m6.c_str());

			myfile >>ue4;
			myfile >>ue5;
			myfile >>ue6;

			ue[0]= atof(ue4.c_str());
			ue[1]= atof(ue5.c_str());
			ue[2]= atof(ue6.c_str());

			myfile >>um4;
			myfile >>um5;
			myfile >>um6;

			um[0]= atof(um4.c_str());
			um[1]= atof(um5.c_str());
			um[2]= atof(um6.c_str());
		
			myfile >> phi45;
			myfile >> phi46;
			myfile >> phi56;

			phi[0]= atof(phi45.c_str());
			phi[1]= atof(phi46.c_str());
			phi[2]= atof(phi56.c_str());

			myfile >>pot;
			myfile >>whicho;
			myfile >>oldchi;
			myfile >>newchi;


			potin= atof(pot.c_str());
			oldchiin= atof(oldchi.c_str());
			newchiin= atof(newchi.c_str());
			whichflagin = atof(whicho.c_str());

	
			neutrinoModel sigModel(mn,ue,um,phi );
			sigModel.numsterile=num_ster;

			for(int i=0;i<(N_e_bins+N_m_bins)*N_dets*vector_modifier; i++){
				myfile >> tempVec;
				vectorin.push_back(atof(tempVec.c_str()));
			}


				if(myfile.eof()){break;}

						
				double chipot = 0;


				//subrract off cosmo
			//	vectorin[0] += -9;
			//	vectorin[N_e_bins+N_m_bins] +=-11;
			//	vectorin[(N_e_bins+N_m_bins)*2] += -10;

				//second position for anti-mode, got to code this better when i have time
				int second = (N_e_bins+N_m_bins)*N_dets;
				
				// subtract off cosmo of anti
			//	if(anti_flag){
			//		vectorin[second] += -9;
			//		vectorin[second+N_e_bins+N_m_bins] += -11;
			//		vectorin[second+(N_e_bins+N_m_bins)*2] += -10;
			//	}

				// now make a mod file, that increases muboone differently
				std::vector<double> mod (N_dets*(N_e_bins+N_m_bins)*vector_modifier, 1.0);
				for(int i =0; i<N_e_bins+N_m_bins; i++){
					mod[i]=ipot;
					mod[i+N_e_bins+N_m_bins]= ipot*0.5+0.5;
					mod[i+2*(N_e_bins+N_m_bins)]=ipot;
				}


				if(anti_flag){
					for(int i =0; i<N_e_bins+N_m_bins; i++){
						mod[i+second]=ipotbar;
						mod[i+second+N_e_bins+N_m_bins] = ipotbar;
						mod[i+second+2*(N_e_bins+N_m_bins)] = ipotbar;
					}
				}


				if(mod.size()!= vectorin.size()){ std::cout<<"WOW big problem"<<std::endl;exit(EXIT_FAILURE);}
				// and modift the prediction
				
				for(int i=0; i<vectorin.size(); i++){
					vectorin[i]=mod[i]*vectorin[i];
		
				}
				//now add back cosmo
			//	vectorin[0]+=9;
			//	vectorin[N_e_bins+N_m_bins]+= 11;
			//	vectorin[(N_e_bins+N_m_bins)*2]+= 10;
				
			//	if(anti_mode == 1){
			//		vectorin[second] +=9;
			//		vectorin[second+N_e_bins+N_m_bins]+=11;
			//		vectorin[second+(N_e_bins+N_m_bins)*2]+=10;
			//	}

			potInstance.calc_chi_POT_vector(sigModel, vectorin, count, ipot, ipotbar); //count and ipot are just there for records

								
			//FAR too much data, dont output this?
			/*
			std::cout<<pred6in[0]*mod[0];
			for(int u=1;u< pred6in.size(); u++){
				std::cout<<" "<<pred6in[u]*mod[u];
			}	
			std::cout<<std::endl;
			*/

			vectorin.clear();
		}//end of file	


	myfile.close();

	}

return 0;
}//end of POT flag



if(inject_flag){

	double ipot =pow(10,pot_num);
	double ipotbar = pow(10,pot_num_bar);

	inputModel.numsterile=num_ster;
	inputModel.phi[0]=inPhi45;


	neutrinoModel testModel = inputModel;

	wrkInstance injectInstance(which_channel, anti_mode, ipot, ipotbar); // anoyingly it has alredy loaded the background model;
	injectInstance.inject_signal(inputModel, which_channel, anti_mode, ipot, ipotbar);


	if(verbose_flag){
		std::cout<<"#v: Inject flag, plotmode : "<<plotmode<<std::endl;
		std::cout<<"#v: Inital POT in nu mode:"<< ipot<<" nu-bar mode: "<<ipotbar<<std::endl;
	         inputModel.printall();
		std::cout<<"#v: #####################################################"<<std::endl;
	}
	if(plotmode == 1)
	{


		injectInstance.init_minim();


		for(double ip = 0; ip <= 2*3.2; ip+=0.1){

			if(margin){
				injectInstance.workingModel = inputModel;
				injectInstance.minimize(ip, ipot, ipotbar);
			}else{
				testModel.phi[0]=ip;
				double ans1=	injectInstance.calc_chi(testModel,1,ipot,ipotbar);
				std::cout<<ip<<" "<<ans1<<std::endl;
			}	
		}// end phi45 loop



	return 0; 
	} // <-- End plotmode==1


	if(plotmode == 2)
	{
		injectInstance.init_minim();

		for(double um4ue4 = -6; um4ue4<0; um4ue4+=0.05){
			for(double ip = 0; ip <= 2*3.2; ip+=0.05){
				if(margin){
					std::cout<<"ERROR: plotmode 2 of inject_signal cannot be run with marginalisation"<<std::endl;
					exit(EXIT_FAILURE);
				}else{
					testModel.phi[0]=ip;
					testModel.Um[0]=1;
					testModel.Ue[0]=sqrt(pow(10,um4ue4));
					testModel.Um[1]=1;
					testModel.Ue[1]=sqrt(pow(10,um4ue4));
					double ans =injectInstance.calc_chi(testModel,1,ipot,ipotbar);
					std::cout<<ip<<" "<<um4ue4<<" "<<ans<<std::endl;
				}	
			}// end phi45 loop
		}// end um4ue4 loop



	return 0; 
	} // <-- End plotmode==2


	if(plotmode ==3)
	{


			//683831 0.870964 1 0 0.045 0.145 0 0.005 0.165 0 5.08938 0 5.08938 1 2 0 69.6192 
		/*
			double Imn[3] = {0.870964,1.0,0};
			double Iue[3] = {0.045,0.115,0};
			double Ium[3] = {0.05,0.125,0};
			double Iphi[3] = {5.0,0.0,0.0};
		
		//.46 .15 .13 .77 .13 .14 5.56
			double Imn[3] = {0.398107,1.0,0};
		//	double Iue[3] = {0.22,0.215,0};
			double Iue[3] = {0.13,0.14,0};
		//	double Ium[3] = {0.24,0.23,0};
			double Ium[3] = {0.15,0.13,0};
			double Iphi[3] = {inPhi45,0.0,0.0};

		
			if(false){ //This is for 1eV^2 signal injection
				Imn[1]=0;
				//Imn[0]=0.660693;// for app
				Imn[0]=1.04713; //for dis
				Iue[1]=0;
				Ium[1]=0;
				Iphi[1]=0;

				//Ium[0]=1;  //for app
				//Iue[0]=0.0570088;	

				Iue[0]=0;
				Ium[0]=0.160182;

			}
		*/

		injectInstance.init_minim();
		neutrinoModel testModel=inputModel;

		//for(double iphi45=0; iphi45<2*3.2; iphi45+=0.1){
			

			double i0 = 0.0;
			double iPi = 3.14159;
			
			if(margin){		
				injectInstance.minimize(i0, ipot, ipotbar);
				injectInstance.minimize(iPi, ipot, ipotbar);
			} else {
				testModel.phi[0]=i0;	
				double ans0 = injectInstance.calc_chi(testModel,1,ipot,ipotbar);
				testModel.phi[0]=iPi;
				double ansPi = injectInstance.calc_chi(testModel,1,ipot,ipotbar);
		
				std::cout<<i0<<" MinimUm: "<<testModel.Ue[0]<<" "<<testModel.Ue[1]<<" "<<testModel.Um[0]<<" "<<testModel.Um[1];
				std::cout<<" "<<inPhi45<<"  " << ans0 << std::endl;
					
				std::cout<<iPi<<" MinimUm: "<<testModel.Ue[0]<<" "<<testModel.Ue[1]<<" "<<testModel.Um[0]<<" "<<testModel.Um[1];
				std::cout<<" "<<inPhi45<<"  " << ansPi << std::endl;
			}
		//} // end of iphi45 True loop THIS loop is too slow for when minimize is on, therfore will do it outside of sbfit.
		return 1;
	}// <-- End Plotmode==3	


	if(plotmode == 4)
	{
	//	injectInstance.init_minim();

			
			for(double ip =0; ip <= 2*3.1501; ip+=0.1){
				
				testModel.phi[0]=ip;
				if(margin){
				//	injectInstance.reset_minim();
				//	injectInstance.init_minim();
				if(verbose_flag)	std::cout<<"justbefore: "<<ip<<" ue4: "<<inputModel.Ue[0]<<" "<<log10(inputModel.Ue[0])<<" ue5: "<<inputModel.Ue[1]<<" "<<log10(inputModel.Ue[1])<<std::endl;
					double ans = injectInstance.minimize(testModel,ipot,ipotbar);
					
					std::cout<<ip<<" "<<pow(inputModel.Ue[0]*inputModel.Um[0],2)<<" "<<pow(inputModel.Ue[1]*inputModel.Um[1],2)<<" "<<ans<<std::endl;
				}else{
					double ans =injectInstance.calc_chi(testModel,1,ipot,ipotbar);
					std::cout<<ip<<" "<<ans<<std::endl;
				}	
			}// end phi45 loop



	return 0; 
	} // <-- End plotmode==2





	// Not sure whats below here. redundant i guess
	//
	int vector_modifier = 1;
	
	if(anti_flag){
		vector_modifier =2;
	}

	/*
	for(double ph=0; ph < 2*3.14145; ph+=0.05){
	double Itest[3]={ph,0,0};

	neutrinoModel testModel(Imn,Iue,Ium,Itest);
		injectInstance.calc_chi(testModel, 0); //count and ipot are just there for records
	}
	return 0;
	*/



/*   ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );
 
   min.SetMaxFunctionCalls(100000);
   min.SetMaxIterations(10000);
   min.SetTolerance(0.01);
 
   ROOT::Math::Functor f(&RosenBrock,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = { -1.,1.2};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   //    min.SetVariable(0,"x",variable[0], step[0]);
   //       min.SetVariable(1,"y",variable[1], step[1]);
   //        
   //           min.Minimize(); 
   //            
   //               const double *xs = min.X();
   //                  cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
   //                          << RosenBrock(xs) << endl;
   //                           
   //                              return 0;
*/





	return 0;
	//now load file with all other vectors.
	std::string filename = "3p1_both";

	if(num_ster == 1){
		switch(which_channel){
			case APP_ONLY:
				filename = "3p1_app";
				break;
			case DIS_ONLY:
				filename = "3p1_dis";
				break;
			case BOTH_ONLY:
				filename = "3p1_both";
				break;
		}

	} else if (num_ster == 2){
			switch(which_channel){
			case APP_ONLY:
				filename = "3p2_app";
				break;
			case DIS_ONLY:
				filename = "3p2_dis";
				break;
			case BOTH_ONLY:
				filename = "3p2_both";
				break;
		}
	} else if(num_ster == 3){
			switch(which_channel){
			case APP_ONLY:
				filename = "3p3_app";
				break;
			case DIS_ONLY:
				filename = "3p3_dis";
				break;
			case BOTH_ONLY:
				filename = "3p3_both";
				break;
		}

	}

	std::string pre =  "fractiondata/";
	std::string post = ".dat";

	if(anti_flag){
		pre = "fractiondata/NUBAR_MODE/";	
		post = ".bar.dat";
	}
	pre.append(filename);
	pre.append(post);
	filename=pre;


if(filename != "none"){

	int k = 0;
	std::string num;
	std::string m4;
	std::string m5;
	std::string m6;

	std::string ue4;
	std::string ue5;
	std::string ue6;

	std::string um4;
	std::string um5;
	std::string um6;

	std::string phi45;
	std::string phi46;
	std::string phi56;

	std::string whicho;

	std::string pot;
	std::string oldchi;
	std::string newchi;

	std::string tempVec;	

	//no longer the pred6in, sometimes pred12 all or whatever bull
	std::vector<double > vectorin;

	double mn[3];
	double ue[3];
	double um[3];
	double phi[3];

	double potin;
	double oldchiin;
	double newchiin;

	int whichflagin = -1;

	int count = 0;

	std::ifstream myfile (filename);
	if (!myfile.is_open())
	{
		std::cout<<"#ERROR: Passed POT file does not exist: "<<filename<<std::endl;
		std::cout<<" which_channel "<<which_channel<<" numster "<<num_ster<<std::endl;
		exit(EXIT_FAILURE);
	}

		while(!myfile.eof()){
		
			count++;	
			myfile >> num;
			//std::cout<<num<<std::endl;	
		
			myfile >> m4;
			myfile >> m5;
			myfile >> m6;
			
			mn[0]= atof(m4.c_str());
			mn[1]= atof(m5.c_str());
			mn[2]= atof(m6.c_str());

			myfile >>ue4;
			myfile >>ue5;
			myfile >>ue6;

			ue[0]= atof(ue4.c_str());
			ue[1]= atof(ue5.c_str());
			ue[2]= atof(ue6.c_str());

			myfile >>um4;
			myfile >>um5;
			myfile >>um6;

			um[0]= atof(um4.c_str());
			um[1]= atof(um5.c_str());
			um[2]= atof(um6.c_str());
		
			myfile >> phi45;
			myfile >> phi46;
			myfile >> phi56;

			phi[0]= atof(phi45.c_str());
			phi[1]= atof(phi46.c_str());
			phi[2]= atof(phi56.c_str());

			myfile >>pot;
			myfile >>whicho;
			myfile >>oldchi;
			myfile >>newchi;


			potin= atof(pot.c_str());
			oldchiin= atof(oldchi.c_str());
			newchiin= atof(newchi.c_str());
			whichflagin = atof(whicho.c_str());

	
			neutrinoModel sigModel(mn,ue,um,phi );
			sigModel.numsterile=num_ster;

			for(int i=0;i<(N_e_bins+N_m_bins)*N_dets*vector_modifier; i++){
				myfile >> tempVec;
				vectorin.push_back(atof(tempVec.c_str()));
			}


				if(myfile.eof()){break;}

						
				double chipot = 0;


				//subrract off cosmo
			//	vectorin[0] += -9;
			//	vectorin[N_e_bins+N_m_bins] +=-11;
			//	vectorin[(N_e_bins+N_m_bins)*2] += -10;

				//second position for anti-mode, got to code this better when i have time
				int second = (N_e_bins+N_m_bins)*N_dets;
				
				// subtract off cosmo of anti
			//	if(anti_flag){
			//		vectorin[second] += -9;
			//		vectorin[second+N_e_bins+N_m_bins] += -11;
			//		vectorin[second+(N_e_bins+N_m_bins)*2] += -10;
			//	}

				// now make a mod file, that increases muboone differently
				std::vector<double> mod (N_dets*(N_e_bins+N_m_bins)*vector_modifier, 1.0);
				for(int i =0; i<N_e_bins+N_m_bins; i++){
					mod[i]=ipot;
					mod[i+N_e_bins+N_m_bins]= ipot*0.5+0.5;
					mod[i+2*(N_e_bins+N_m_bins)]=ipot;
				}


				if(anti_flag){
					for(int i =0; i<N_e_bins+N_m_bins; i++){
						mod[i+second]=ipotbar;
						mod[i+second+N_e_bins+N_m_bins] = ipotbar;
						mod[i+second+2*(N_e_bins+N_m_bins)]=ipotbar;
					}
				}


				if(mod.size()!= vectorin.size()){ std::cout<<"WOW big problem"<<std::endl;exit(EXIT_FAILURE);}
				// and modift the prediction
				
				for(int i=0; i<vectorin.size(); i++){
					vectorin[i]=mod[i]*vectorin[i];
		
				}
				//now add back cosmo
			//	vectorin[0]+=9;
			//	vectorin[N_e_bins+N_m_bins]+= 11;
			//	vectorin[(N_e_bins+N_m_bins)*2]+= 10;
				
			//	if(anti_mode == 1){
			//		vectorin[second] +=9;
			//		vectorin[second+N_e_bins+N_m_bins]+=11;
			//		vectorin[second+(N_e_bins+N_m_bins)*2]+=10;
			//	}

			injectInstance.calc_chi_POT_vector(sigModel, vectorin, count, ipot, ipotbar); //count and ipot are just there for records

								
			//FAR too much data, dont output this?
			/*
			std::cout<<pred6in[0]*mod[0];
			for(int u=1;u< pred6in.size(); u++){
				std::cout<<" "<<pred6in[u]*mod[u];
			}	
			std::cout<<std::endl;
			*/

			vectorin.clear();
		}//end of file	


	myfile.close();

	}


}//end of injection flag






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
	std::vector<std::vector<double >> vMcI = to_vector(McI);
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
//		samplespectrum.SetNuBarMode();

		samplespectrum.oscillate_sample();
		break;
	}

	

}

if(cov_flag){


}
 
if(sens_flag){
	std::cout<<"Starting 3p1 sensitivity scan"<<std::endl;
        wrkInstance sensInstance(which_channel, anti_mode);

	for(double m = -2.00; m <=2.04; m=m+0.04){
		

			for(double sins2 = log10(0.25) ; sins2 > -5; sins2 = sins2 - 0.2){
				double t_ue = -99;
				double t_um = -99;
				double t_chi = 1e5;


			for(double uei2 = log10(0.3) ; uei2 > -10; uei2 = uei2 - 0.05){
				double uei = pow(10,uei2);
				
				double umi = sqrt(pow(10,sins2))/(2*uei);
				if(umi > 0.3){continue;};
			
				double imn[3] = {sqrt(pow(10,m)),0,0};
				double iue[3] = {umi,0,0};
				double ium[3] = {uei,0,0};
				double iph[3] = {0,0,0};

				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile = 1;
				sensInstance.isVerbose = false;
	

				double chi2A = sensInstance.calc_chi(signalModel, 0);
				if(chi2A < t_chi){
					t_chi =chi2A;
					t_ue = umi;
					t_um = uei;
			//		if(which_channel==APP_ONLY){break;}	
				}


				 iue[0]= uei;
				 ium[0] = umi;

				neutrinoModel signalModel2(imn,iue,ium,iph);
				signalModel2.numsterile = 1;
				
				chi2A = sensInstance.calc_chi(signalModel2, 0);
				if(chi2A < t_chi){
					t_chi =chi2A;
					t_ue = uei;
					t_um = umi;
			//		if(which_channel==APP_ONLY){break;}	
				}
				
	//	std::cout<<m<<" "<<sins2<<" "<<chi2A<<" "<<uei<<" "<<umi<<std::endl;		
	}//end ue

		std::cout<<m<<" "<<sins2<<" "<<t_chi<<" "<<t_ue<<" "<<t_um<<std::endl;		
	}//end sins2
	}//end m





return 0;
}


// OLD sensitivities...
if(sens_flag && false){

	double norm_pot = 1.0;
	wrkInstance signalInstance(which_channel , anti_mode);


				for(double phi45=0; phi45<2*3.14159; phi45+=0.1){
	
					double imn[3] = {1,10,0};
					double iue[3] = {0.1,0.1,0};
					double ium[3] = {0.1,0.1, 0};
					double iph[3] ={phi45,0,1.0};

					neutrinoModel signalModel(imn,iue,ium,iph);
					signalModel.numsterile= 2;

					signalInstance.calc_chi(signalModel, 0);

				}

}


if(sens_flag && false)
{

 	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	SBN_detector * SBND_mu = new SBN_detector(0,true);
 	SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;

	if(dis_flag && !app_flag){
		//usedetsys=false;
	}

	neutrinoModel nullModel;
	SBN_spectrum bkgspec(nullModel);
	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
			if(false){
				double modd= 0.5;	
				for(int i = 0; i < N_e_bins; i++){

				double pot =1e-12;	
					bkgspec.sbnd_e[i]= bkgspec.sbnd_e[i]*pot;
					bkgspec.sbnd_e_pho[i]= bkgspec.sbnd_e_pho[i]*pot;
					bkgspec.sbnd_e_dirt[i]= bkgspec.sbnd_e_dirt[i]*pot;
					bkgspec.sbnd_e_mu[i]= bkgspec.sbnd_e_mu[i]*pot;
					bkgspec.sbnd_e_cosmo[i]=pot;
					bkgspec.sbnd_f[i]= bkgspec.sbnd_f[i]*pot;
					bkgspec.sbnd_f_bar[i]= bkgspec.sbnd_f_bar[i]*pot;
		
					bkgspec.icarus_e[i]= bkgspec.icarus_e[i]*pot;
					bkgspec.icarus_e_pho[i]= bkgspec.icarus_e_pho[i]*pot;
					bkgspec.icarus_e_dirt[i]= bkgspec.icarus_e_dirt[i]*pot;
					bkgspec.icarus_e_mu[i]= bkgspec.icarus_e_mu[i]*pot;
					bkgspec.icarus_e_cosmo[i]=pot;
					bkgspec.icarus_f[i]= bkgspec.icarus_f[i]*pot;
					bkgspec.icarus_f_bar[i]= bkgspec.icarus_f_bar[i]*pot;
							
					bkgspec.uboone_e[i]= bkgspec.uboone_e[i]*modd;
					bkgspec.uboone_e_pho[i]= bkgspec.uboone_e_pho[i]*modd;
					bkgspec.uboone_e_dirt[i]= bkgspec.uboone_e_dirt[i]*modd;
					bkgspec.uboone_e_mu[i]= bkgspec.uboone_e_mu[i]*modd;
					bkgspec.uboone_f[i]= bkgspec.uboone_f[i]*modd;
					bkgspec.uboone_f_bar[i]= bkgspec.uboone_f_bar[i]*modd;

				}
			
				for(int i =0; i< N_m_bins; i++){
					double pot = 1e-12;
					bkgspec.sbnd_m[i]= bkgspec.sbnd_m[i]*pot;
					bkgspec.sbnd_m_pion[i]= bkgspec.sbnd_m_pion[i]*pot;

					bkgspec.uboone_m[i]= bkgspec.uboone_m[i]*modd;
					bkgspec.uboone_m_pion[i]= bkgspec.uboone_m_pion[i]*modd;

					bkgspec.icarus_m[i]= bkgspec.icarus_m[i]*pot;
					bkgspec.icarus_m_pion[i]= bkgspec.icarus_m_pion[i]*pot;

				}
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
	
				AppSpec.load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
				AppSpec.load_freq_3p3(SBND);
				AppSpec.load_freq_3p3(UBOONE);


					

				for(int i = 0; i < N_e_bins; i++){
				double modd = 0.5;
				double pot =1e-12;	
					AppSpec.sbnd_e[i]= AppSpec.sbnd_e[i]*pot;
					AppSpec.sbnd_e_pho[i]= AppSpec.sbnd_e_pho[i]*pot;
					AppSpec.sbnd_e_dirt[i]= AppSpec.sbnd_e_dirt[i]*pot;
					AppSpec.sbnd_e_mu[i]= AppSpec.sbnd_e_mu[i]*pot;
					AppSpec.sbnd_e_cosmo[i]=pot;
					AppSpec.sbnd_f[i]= AppSpec.sbnd_f[i]*pot;
					AppSpec.sbnd_f_bar[i]= AppSpec.sbnd_f_bar[i]*pot;
		
					AppSpec.icarus_e[i]= AppSpec.icarus_e[i]*pot;
					AppSpec.icarus_e_pho[i]= AppSpec.icarus_e_pho[i]*pot;
					AppSpec.icarus_e_dirt[i]= AppSpec.icarus_e_dirt[i]*pot;
					AppSpec.icarus_e_mu[i]= AppSpec.icarus_e_mu[i]*pot;
					AppSpec.icarus_e_cosmo[i]=pot;
					AppSpec.icarus_f[i]= AppSpec.icarus_f[i]*pot;
					AppSpec.icarus_f_bar[i]= AppSpec.icarus_f_bar[i]*pot;

	
					AppSpec.uboone_e[i]= AppSpec.uboone_e[i]*modd;
					AppSpec.uboone_e_pho[i]= AppSpec.uboone_e_pho[i]*modd;
					AppSpec.uboone_e_dirt[i]= AppSpec.uboone_e_dirt[i]*modd;
					AppSpec.uboone_e_mu[i]= AppSpec.uboone_e_mu[i]*modd;
					AppSpec.uboone_f[i]= AppSpec.uboone_f[i]*modd;
					AppSpec.uboone_f_bar[i]= AppSpec.uboone_f_bar[i]*modd;

								
				}
			
				for(int i =0; i< N_m_bins; i++){
					double pot = 1e-12;
				double modd = 0.5;
					AppSpec.sbnd_m[i]= AppSpec.sbnd_m[i]*pot;
					AppSpec.sbnd_m_pion[i]= AppSpec.sbnd_m_pion[i]*pot;

					AppSpec.uboone_m[i]= AppSpec.uboone_m[i]*modd;
					AppSpec.uboone_m_pion[i]= AppSpec.uboone_m_pion[i]*modd;

					AppSpec.icarus_m[i]= AppSpec.icarus_m[i]*pot;
					AppSpec.icarus_m_pion[i]= AppSpec.icarus_m_pion[i]*pot;

				}


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
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}


				double sin22em = 4.0*pow(appearanceModel.Ue[0]*appearanceModel.Um[0],2.0);

				std::cout<<m<<" "<<appearanceModel.Ue[0]<<" "<<appearanceModel.Um[0]<<" "<<chi2<<" "<<sin22em<<std::endl;
			}//end random u run
		}//end mass run
	} //end 3p1 APPearance only sensitivity analysis

	if(sens_num == 1 && dis_flag && dis_which == 1)
	{
		std::cout<<"Begining  N=1, muon dis only"<<std::endl;
		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(int i = 0; i< 300; i++){

				double umi = rangen->Uniform(-0.14,-2.0);
				neutrinoModel disappearanceModel(sqrt(pow(10,m)), 0.0, pow(10,umi));

				disappearanceModel.dm41Sq = pow(10,m);

				SBN_spectrum DisSpec(disappearanceModel);
				
				DisSpec.load_freq_3p3(UBOONE);				
				DisSpec.load_freq_3p3(ICARUS);
				DisSpec.load_freq_3p3(SBND);
				
	
			//	DisSpec.load_bkg(ICARUS);
			//	DisSpec.load_bkg(SBND);


				for(int i = 0; i < N_e_bins; i++){
				double modd = 0.5;
				double pot =1e-12;	
					DisSpec.sbnd_e[i]= DisSpec.sbnd_e[i]*pot;
					DisSpec.sbnd_e_pho[i]= DisSpec.sbnd_e_pho[i]*pot;
					DisSpec.sbnd_e_dirt[i]= DisSpec.sbnd_e_dirt[i]*pot;
					DisSpec.sbnd_e_mu[i]= DisSpec.sbnd_e_mu[i]*pot;
					DisSpec.sbnd_e_cosmo[i]=pot;
					DisSpec.sbnd_f[i]= DisSpec.sbnd_f[i]*pot;
					DisSpec.sbnd_f_bar[i]= DisSpec.sbnd_f_bar[i]*pot;
		
					DisSpec.icarus_e[i]= DisSpec.icarus_e[i]*pot;
					DisSpec.icarus_e_pho[i]= DisSpec.icarus_e_pho[i]*pot;
					DisSpec.icarus_e_dirt[i]= DisSpec.icarus_e_dirt[i]*pot;
					DisSpec.icarus_e_mu[i]= DisSpec.icarus_e_mu[i]*pot;
					DisSpec.icarus_e_cosmo[i]=pot;
					DisSpec.icarus_f[i]= DisSpec.icarus_f[i]*pot;
					DisSpec.icarus_f_bar[i]= DisSpec.icarus_f_bar[i]*pot;

	
					DisSpec.uboone_e[i]= DisSpec.uboone_e[i]*modd;
					DisSpec.uboone_e_pho[i]= DisSpec.uboone_e_pho[i]*modd;
					DisSpec.uboone_e_dirt[i]= DisSpec.uboone_e_dirt[i]*modd;
					DisSpec.uboone_e_mu[i]= DisSpec.uboone_e_mu[i]*modd;
					DisSpec.uboone_f[i]= DisSpec.uboone_f[i]*modd;
					DisSpec.uboone_f_bar[i]= DisSpec.uboone_f_bar[i]*modd;

								
				}
			
				for(int i =0; i< N_m_bins; i++){
					double pot = 1e-12;
				double modd = 0.5;
					DisSpec.sbnd_m[i]= DisSpec.sbnd_m[i]*pot;
					DisSpec.sbnd_m_pion[i]= DisSpec.sbnd_m_pion[i]*pot;

					DisSpec.uboone_m[i]= DisSpec.uboone_m[i]*modd;
					DisSpec.uboone_m_pion[i]= DisSpec.uboone_m_pion[i]*modd;

					DisSpec.icarus_m[i]= DisSpec.icarus_m[i]*pot;
					DisSpec.icarus_m_pion[i]= DisSpec.icarus_m_pion[i]*pot;

				}


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


				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += (back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
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
						chi2 += mod*(back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}


				double sin22ee = 4.0*(1-pow(disappearanceModel.Ue[0],2.0))*pow(disappearanceModel.Ue[0],2.0);

				std::cout<<m<<" "<<disappearanceModel.Ue[0]<<" "<<chi2<<" "<<sin22ee<<std::endl;
			}//end random ue4
		}//end m for loop
	} //end 3p1 sensitivity disapearance (ue4 dis only) only analysis 




	if(sens_num == 1 && both_flag && mode_flag == 1)
	{
	std::cout<<"Begining N=1, dual appearance and dissapearance, optimised for sin^2them   : 3p1_both_em.dat"<<std::endl;

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
     

		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(double sins2 = log10(0.25) ; sins2 > -5; sins2 = sins2 - 0.2){
			for(double uei2 = log10(0.4) ; uei2 > -3; uei2 = uei2 - 0.05){
				double uei = pow(10,uei2);
				
				double umi = sqrt(pow(10,sins2))/(2*uei);
				if(umi > 0.4){continue;};

				neutrinoModel bothModel(sqrt(pow(10,m)), uei , umi);
				bothModel.dm41Sq = pow(10,m);
			
				
			//	double imn[3] = {sqrt(pow(10,m)),sqrt(pow(10,1.24)),0.0};
			//	double iue[3] = {uei,0.069,0};
			//	double ium[3] = {umi, 0.16, 0.0};
			//	double iph[3] = {1.8*3.14159,0.0, 0.0};

			//	neutrinoModel bothModel(imn,iue,ium,iph);
				


				SBN_spectrum BothSpec(bothModel);
				BothSpec.which_mode = 2;

				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

	

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

				ntuple.Fill(bothModel.dm41Sq,bothModel.Ue[0],bothModel.Um[0],chi2);
				std::cout<<bothModel.dm41Sq<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<pow(10,sins2)<<std::endl;
			}//end um4
			std::cout<<"#Finished m: "<<m<<" "<<sins2<<std::endl;
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

if(sens_num == 2)   // Fix everything global, vary phi
	{
	std::cout<<"Begining N=2, fix all vary phi "<<std::endl;



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
     

			//for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.075){
			//for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.075){
			

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
			
				double minPhi = 0;
				double minChi = 10000;	
				for(double iphi = 0; iphi< 2; iphi+=0.01){
				//double iphi = 1.8;
			
				
				//0.69, 1.3
				double imn[3] = {sqrt(pow(10,-0.16)),sqrt(pow(10,0.12)),0.0};
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
				std::cout<<chi2<<" "<<iphi*3.14159<<std::endl;
				//std::cout<<m4<<" "<<m5<<" "<<bothModel.dm54Sq<<" "<<std::endl;
				} // end phi 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis
	
if(sens_num == 2&& false )   // This is the m41 m51 margined phi case for best fit
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
		}//en  for loop

 
		
//		outputFile.cd();
//		std::cout<<"write ntuple"<<std::endl;
//		ntuple.Write();
//	 	std::cout<<"close file"<<std::endl;
  // 		outputFile.Close();
//		std::cout<<"end all"<<std::endl;
	} //end 3p2 sensitivity both analysis
}//end sens_flag

	


if(anti_flag && false){ //depriciated, now include in wrkInstance

	if(verbose_flag) std::cout<<"#**********Begining nu+nubar mode analysis**********"<<std::endl;
	SBN_detector * ICARUS = new SBN_detector(2);
 	SBN_detector * SBND = new SBN_detector(0);
 	SBN_detector * UBOONE = new SBN_detector(1);

	//SBN_detector * ICARUS_mu = new SBN_detector(2,true);
 	//SBN_detector * SBND_mu = new SBN_detector(0,true);
 	//SBN_detector * UBOONE_mu = new SBN_detector(1,true);
	bool usedetsys = true;

	if(dis_flag && !app_flag){
		usedetsys=false;
	}

	neutrinoModel nullModel;


	if(verbose_flag) std::cout<<"# Initialising nu+nubar mode backgrrounds"<<std::endl;
	
	SBN_spectrum bkgspec(nullModel);
	SBN_spectrum bkgbarspec(nullModel);
	bkgbarspec.SetNuBarMode();

	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);
			
	bkgbarspec.load_bkg(ICARUS);
	bkgbarspec.load_bkg(SBND);
	bkgbarspec.load_bkg(UBOONE);
	
	std::vector<double > back6 = bkgspec.get_sixvector();
	std::vector<double > back9 = bkgspec.get_ninevector();
	std::vector<double > back  = bkgspec.get_vector();

	std::vector<double > backbar6 = bkgbarspec.get_sixvector();
	std::vector<double > backbar9 = bkgbarspec.get_ninevector();
	std::vector<double > backbar  = bkgbarspec.get_vector();

	std::vector<double> back_all = back;
	back_all.insert(back_all.end(), backbar.begin(), backbar.end() );

	std::vector<double> back_all_12 =back6;
       	back_all_12.insert(back_all_12.end(), backbar6.begin(), backbar6.end() );


	TRandom *rangen    = new TRandom();


		int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti;
		int contMsize = (N_e_bins+N_m_bins)*N_dets*N_anti;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */

		TMatrixT <double> McI(contMsize, contMsize);
		TMatrixT <double> McIbar(contMsize, contMsize);

		if(verbose_flag) std::cout<<"# Filling up Systematics correlation matrix, size: "<<bigMsize <<std::endl;
		// Fill systematics from pre-computed files
		TMatrixT <double> Msys(bigMsize,bigMsize);
		sys_fill(Msys,usedetsys);

		// systematics per scaled event
		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back_all[i]*back_all[j];
			}
		}




		if(verbose_flag) std::cout<<"# Filling up Statistics correlation matrix, size: "<<bigMsize <<std::endl;
		// Fill stats from the back ground vector
		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back_all);

		//And then define the total covariance matrix in all its glory
		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}

		// Now contract back the larger antimatrix
		TMatrixT<double > Mctotal(contMsize,contMsize);

		if(verbose_flag) std::cout<<"# Begin nu+nubar mode matrix contraction: contract_signal2()"<<std::endl;
		contract_signal2_anti(Mtotal,Mctotal);
	
		if(verbose_flag) std::cout<<"# Finished contraction, matrix size reduced from: "<<Mtotal.GetNcols()<<" to: "<<Mctotal.GetNcols()<<std::endl;

		// just to hold determinant
		double invdet=0; 

		// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
		McI = Mctotal.Invert(&invdet);


		// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
		std::vector<std::vector<double >> vMcI = to_vector(McI);


				double m4 = 0.04;
				double m5 = -0.04;			

				//neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				//bothModel.dm41Sq = pow(10,m);
			
				double iphi = 1.58;	
				
				double imn[3] = {sqrt(pow(10,m4)),sqrt(pow(10,m5)),0.0};
				double iue[3] = {0.15,0.13,0};  // These are best fit!
				double ium[3] = {0.069,0.16, 0.0};
				//double iue[3] = {0.2,0.2,0};  // generic ones
				//double ium[3] = {0.2,0.2,0.0};


				double iph[3] = {iphi*3.14159, 0.0, 0.0};
				
				
				neutrinoModel bothModel(imn,iue,ium,iph);
				bothModel.numsterile =2;

				double round54 = round(log10(fabs(bothModel.dm54Sq))/0.04)*0.04;

			//	if(fabs(bothModel.dm54Sq) >= 100){ 
				//	std::cout<<"skipping this one 1:"<<std::endl;
			//			continue;
			//	}
			//	if(round54 > 2 ){ 
			//		std::cout<<"skipping this one 1: round54 "<<round54<<std::endl;
			//			continue;
			//	}

//				std::cout<<"dm54: "<<bothModel.dm54Sq<<" "<<bothModel.dm64Sq<<" "<<bothModel.dm65Sq<<std::endl;
				SBN_spectrum BothSpec(bothModel);
				SBN_spectrum BothBarSpec(bothModel);
				BothBarSpec.SetNuBarMode();

				BothSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND);
				BothSpec.load_freq_3p3(UBOONE);

				BothBarSpec.load_freq_3p3(ICARUS);//1 is stupid dis flag (temp)
				BothBarSpec.load_freq_3p3(SBND);
				BothBarSpec.load_freq_3p3(UBOONE);
			
		
		
				std::vector<double > pred6 = BothSpec.get_sixvector();
				std::vector<double > predbar6 = BothBarSpec.get_sixvector();
		
				std::vector<double> pred_all_12 = pred6;	
			     	pred_all_12.insert( pred_all_12.end(), predbar6.begin(), predbar6.end() );
		
				double chi2=0;
		
			
				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

				int whatsize = McI.GetNcols();



				if( McI.GetNcols() != McI.GetNrows() && back_all_12.size() != pred_all_12.size() && back_all_12.size() != McI.GetNcols() )
				{
					std::cout<<"ERROR: Something is wrong. In Chi^2 calc matrix != vector length "<<std::endl;
				}



				if(verbose_flag) std::cout<<"# calculate chi^2 value using contracted matrix!"<<std::endl;

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						chi2 += (back_all_12[i]-pred_all_12[i])*vMcI[i][j]*(back_all_12[j]-pred_all_12[j]);
					}
				}
			
			

				//std::cout<<m<<" "<<bothModel.Ue[0]<<" "<<bothModel.Um[0]<<" "<<chi2<<" "<<std::endl;
				/*nCHI = chi2;
				nUM4 = bothModel.Um[0];
				nUE4 = bothModel.Ue[0];
				nDM4 = bothModel.dm41Sq;*/
//				ntuple.Fill(pow(10,m4),pow(10,m5),pow(10,ueiMin),pow(10,umiMin),chiMin);
				std::cout<<pow(10,m4)<<" "<<pow(10,m5)<<" "<<chi2<<std::endl;


}






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

	bkgspec.neutral_test(SBND);
//	bkgspec.neutral_test(SBND);
//	bkgspec.neutral_test(SBND);


return 0;	
	bkgspec.load_bkg(ICARUS);
	bkgspec.load_bkg(SBND);
	bkgspec.load_bkg(UBOONE);


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




