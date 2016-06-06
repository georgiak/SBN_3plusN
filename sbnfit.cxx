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

int  sens_num=1;
double dm41 = -1.0;

double in_dm = 0;
double in_ue4 = 0;
double in_um4=0;

bool stat_only = false;
int dis_which = 1;

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
	{"stat-only",		no_argument,		0, 'f'},
	{"sample",		no_argument,		0, 's'},
	{"cov",			no_argument, 		0, 'c'},
	{"dis",			required_argument, 	0, 'd'},
	{"app",			no_argument,		0, 'a'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "d:afu:e:m:svhS:cFTB", longopts, &index);

	switch(iarg)
	{
		case 'F':
			fit_flag = true;
			//mS = strtof(optarg,NULL);
			break;
		case 'B':
			bkg_flag = true;
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
		case 'f':
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
			std::cout<<"\t-d\t--app\t\tRun app only sensitivity"<<std::endl;
			std::cout<<"\t-a\t--dis\t\tRun dis only sensitivity, takes 1 arg 0 for e-dis 1 for mudis"<<std::endl;
			std::cout<<"\t--stat-only\t\t Run with no systematics, ony stats"<<std::endl;
			std::cout<<"\t-B\t--bkg\t\trun bkg generating test code"<<std::endl;
			std::cout<<"\t-h\t--help\t\tDisplays this help message"<<std::endl;
			std::cout<<"\t-s\t--sample\t\t generate all the sin and sin^2 samples"<<std::endl;
			std::cout<<"\t-d\t\t\tRequired Argument. Creates a sin and sin^2 frequency ntuples for a dmsq."<<std::endl;
			std::cout<<"\t-v\t\t\tVerbose run, mostly debugging"<<std::endl;	
			return 0;
	}

}

if(app_flag&&dis_flag){both_flag = true; app_flag=false; dis_flag=false;}


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

	system("mv bkg_data/NUBAR_MODE/ICARUS_-inf.root bkg_data/NUBAR_MODE/ICARUS_bkg.root"); 
	system("mv bkg_data/NUBAR_MODE/SBND_-inf.root bkg_data/NUBAR_MODE/SBND_bkg.root"); 
	system("mv bkg_data/NUBAR_MODE/uBooNE_-inf.root bkg_data/NUBAR_MODE/uBooNE_bkg.root");

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
	
	for(double m = -2.0; m <=2; m=m+0.04){
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
     

		for(double m = -2.00; m <=2.04; m=m+0.04){
			for(double umi = log10(0.5); umi >= -3.0; umi = umi - 0.05){
			for(double uei = log10(0.5); uei >= -3.0; uei = uei - 0.05){
				
				neutrinoModel bothModel(sqrt(pow(10,m)), pow(10,uei),pow(10,umi));
				bothModel.dm41Sq = pow(10,m);


				SBN_spectrum BothSpec(bothModel);
				
				BothSpec.load_freq_3p3(ICARUS_mu);//1 is stupid dis flag (temp)
				BothSpec.load_freq_3p3(SBND_mu);
				BothSpec.load_freq_3p3(UBOONE_mu);

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





