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
 ************************************************************
 ************************************************************/

#define N_m_bins 20
#define N_e_bins 11
#define N_dets 3

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


double dm41 = -1.0;
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
	{"help",		no_argument, 		0, 'h'},
	{"verbose",		no_argument, 		0, 'v'},
	{"dm",			required_argument, 	0, 'd'},
	{0,			no_argument, 		0,  0},
};

int index; 
int iarg = 0;
opterr=1;

while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "d:vhFTB", longopts, &index);

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
		case 'd':
			dm41 = strtof(optarg,NULL); 
			break;
		case '?':
		case 'h':
			std::cout<<"Allowed arguments:"<<std::endl;
			std::cout<<"\t-F\t--fit\t\trun SBN fitting code"<<std::endl;
			std::cout<<"\t-T\t--test\t\trun SBN test code"<<std::endl;
			std::cout<<"\t-B\t--bkg\t\trun bkg generating test code"<<std::endl;
			std::cout<<"\t-h\t--help\t\tDisplays this help message"<<std::endl;
			std::cout<<"\t-d\t\t\tRequired Argument. run a 3+1 with dm41."<<std::endl;
			std::cout<<"\t-v\t\t\tVerbose run, mostly debugging"<<std::endl;	
			return 0;
	}

}


//Begin program flow control
if(fit_flag){

/*************************************************************
 *************************************************************
 *		Begin SBN Fitting Code
 ************************************************************
 ************************************************************/

	//initialize random generator	
	TRandom *rangen    = new TRandom();

	// Calculate the expected background. Will read this in once trusted. see bkg_flag code below for how this works.
	neutrinoModel nullModel;
	SBN_spectrum bkgspectrum(nullModel);
	bkgspectrum.oscillate();

	std::vector<double > back6 = bkgspectrum.get_sixvector();
	std::vector<double > back9 = bkgspectrum.get_ninevector();

	/*Create a 3+1 spectrum object
	 */
	double m4 = 2.0; //in eV
	double ue4 = 0.01;
	double um4 = 0.02;

	neutrinoModel sterilemodel(m4, ue4, um4);
	SBN_spectrum SBN3p1(sterilemodel);

	/*Can do basically all I need using this structure
	 *  see model.h for more details
	 *
	 * For now lets calculate the expected spectra of the whole SBN network
	 * nu_mu -> nu_mu disapearance at SBND, uBooNE and ICARUS
	 * nu_e intrinsic -> nu_e intrinisc disapearance at ''
	 * nu_mu -> nu_e appearance at ''
	 *
	 * Does NOT calculate nu_e intrinsic -> nu_mu appearace (although could easily actually)
	 *
	 * These can be called individually by
	 *	sterilemodel.fill_app(SBN_detector * SBND);
	 *	sterilemodel.fill_intrin(SBN_detector * ICARuS);
	 *	.. etc..
	 * However they are bundeled together in the SBN_spectrum::oscillate() functionality
	 *  */

	SBN3p1.oscillate(); 

	std::vector<double > pred6 = SBN3p1.get_sixvector();
	std::vector<double > pred9 = SBN3p1.get_ninevector();

		/* DEPRICIATED IGNORE
		double s2 = rangen->Uniform(-5,-2);
		std::cout<<"# seed "<<s2<<std::endl;
		double max = log10(sqrt(pow(10.0,s2))/2.0);
		double um = pow(10.0, rangen->Uniform(max,0.0));
		double ue = sqrt(pow(10.0,s2))/(2.0*um);
		double Sapp=4*pow(testModel.Ue[0]*testModel.Um[0],2.0);
		double Sdis=1.0-4*pow(testModel.Um[0],2)*(1- pow(testModel.Um[0],2));*/


	/******************************************************
	 *	Actually calculating the fit will eventually be 
	 *	its own class soon. currently ad hoc.
	 *
	 * ****************************************************/

	int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
	int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

	/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
	 * */
	TMatrixT <double> M(matrix_size,matrix_size);
	TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
	TMatrixT <double> McI(matrix_size_c, matrix_size_c);

	//Just assume stats only at this point.
	stats_fill(M,back9);
	//M.Print();
	
	double invdet=0; // just to hold determinant
	//MI=M.Invert(&invdet);

	contract_signal(M,Mc);
	
	//Mc.Print();
	
	double chi2=0;
		
	//	bit o inverting, root tmatrix seems perfectly fast	
	McI = Mc.Invert(&invdet);

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

	std::cout<<" m4: "<<sterilemodel.mNu[0]<<" Ue4:  "<<sterilemodel.Ue[0] <<" Um4: "<<sterilemodel.Um[0]<<" chi2: "<<chi2<<std::endl;

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
	neutrinoModel nullModel;

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
	bkgspectrum.oscillate();

	//and keep the end histograms in internal variable, sbnd_e sbnd_m ..etc..	
	//as well as writing them all to root files for plotting.. etc..
	// NO plotting should be done in this sbnfit code. doesnt seem right. 

	std::cout<<"****************** mu-spectra *******************************"<<std::endl;
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


if(test_flag){
	/*************************************************************
	 *************************************************************
	 *		test flag please ignore, 
	 *	no really, just a jumble place for running tests
	 ************************************************************
	 ************************************************************/

 	SBN_detector * ICARUS = new SBN_detector(2);

	double mn[3] = {0.8,0.2,12};
	double ue[3] = {0.5,0.13,0.2};
	double um[3] = {0.5,0.3,0.2};
	double ph[3] = {0.0,0.0,0.0};

	neutrinoModel testModel(mn,ue,um,ph);
	neutrinoModel testModel2(2.0, ue[0], um[0]);

//	for(double ee =0.01; ee< 5; ee=ee+0.005){
//		std::cout<<ee<<" "<<testModel.oscProb(2,1,ee,0.6)<<std::endl;
//	}

	SBN_spectrum myspec(testModel2);

	double Sapp=4*pow(testModel.Ue[0]*testModel.Um[0],2.0);
	double Sdis=1.0-4*pow(testModel.Um[0],2)*(1- pow(testModel.Um[0],2));
	std::cout<<"Ptest2 "<<Pmue(600,1.0,1.0,Sdis)<<std::endl;


	myspec.oscillate();
	myspec.vec_print();







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





