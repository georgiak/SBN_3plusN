#ifndef MODEL_H_
#define MODEL_H_

#include <cmath>
#include <vector>
#include "TH1F.h"
#include <iostream>
#include "TRandom.h"
#include <algorithm>
#include "detector.h"
#include "prob.h"

/*********************************************
* Data Struct: Same as Davio's output.
* -------------------------------------------
* ToDo: 
* Should be packed in a nTuple to me later. 
*
*Get the intrinsic nu_e to oscillate toooo
*
* *******************************************/
#define MPROTON  0.938
#define MPION   0.13957
#define MKAON   0.493
#define MSIGMA  1.189

#define psmear  0.05
#define pismear  0.05
#define EMsmear  0.15
#define MUsmear  0.06

#define p_thresh  0.021
#define pip_thresh  0.021
#define pim_thresh  0.021
#define vertex_thresh  0.05
#define EM_thresh  0.200



/*********************************************
* Main Class  SBN_spectrum
* -------------------------------------------
*	main workhorse of the code
*
* ToDo:
*	class needs appearance and disapearance fiducial volumes! they arent the same  
* 
* *******************************************/

class SBN_spectrum {
	struct neutrinoModel nullModel;
	static const int N_e_bins = 11;
	static const int N_m_bins = 19;
	static constexpr double mu_bins[N_m_bins+1]=   {0.2,0.3,0.4,0.45,0.5,0.55,0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.25, 1.5, 2., 2.5, 3.};
	static constexpr double e_bins[N_e_bins+1]= {0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.3,1.5,1.75,2,3};	 
	static constexpr double ebinw[N_e_bins]   = {0.15,0.15,0.15,0.15,0.15,0.15,0.2,0.2,0.25,0.25,1.0};
	static constexpr double mubinw[N_m_bins] = {0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.25, 0.25, 0.5, 0.5, 0.5};

	public:
	std::vector<double >	uboone_e ;
	std::vector<double > 	sbnd_e ;
	std::vector<double >    icarus_e ;

	std::vector<double >  	sbnd_m ;
  	std::vector<double >  	uboone_m ;
	std::vector<double >  	icarus_m  ;

	std::vector<double >  	sbnd_m_pion ;
  	std::vector<double >  	uboone_m_pion ;
	std::vector<double >  	icarus_m_pion  ;

	std::vector<double >	uboone_f ;
	std::vector<double > 	sbnd_f ;
	std::vector<double >    icarus_f ;

	std::vector<double > sbnd_e_pho;
	std::vector<double > uboone_e_pho;
	std::vector<double > icarus_e_pho;

	std::vector<double > sbnd_e_mu;
	std::vector<double > uboone_e_mu;
	std::vector<double > icarus_e_mu;

	std::vector<double > sbnd_e_dirt;
	std::vector<double > uboone_e_dirt;
	std::vector<double > icarus_e_dirt;

	std::vector<double > sbnd_e_cosmo;
	std::vector<double > uboone_e_cosmo;
	std::vector<double > icarus_e_cosmo;

	std::vector<double > sbnd_f_bar;
	std::vector<double > uboone_f_bar;
	std::vector<double > icarus_f_bar;

	//Redundant histograms here really.
	TH1D * THuboone_e,* THuboone_m, * THsbnd_e, * THsbnd_m,* THicarus_e,* THicarus_m;

	struct neutrinoModel workingModel;	
	
	//Constructors, if blank corresponds to background only. 
	SBN_spectrum (struct neutrinoModel nuModel);
	SBN_spectrum ();

	std::vector<double > get_ninevector();
	std::vector<double > get_sixvector();
	
	void update_model(struct neutrinoModel nuModel);
	
	std::vector<double > get_vector();
	

	int fill_app(SBN_detector *);		//Somewhat depriciated functions
	int fill_dis(SBN_detector *);
	int fill_intrin(SBN_detector *);


	int fill_app_sample(SBN_detector *); 	//these are really the ones that are used 
	int fill_dis_sample(SBN_detector *);
	int fill_intrin_sample(SBN_detector *);

	std::vector<double > add_SBN_spectrum(SBN_spectrum other);

	void oscillate();
	void oscillate_sample();
	
	void fill_hists();
        void fill_vectors();	
	void vec_print();

	int load_freq(SBN_detector *, int);
	int load_freq_3p3(SBN_detector *);
 	int load_bkg(SBN_detector *); 


	double prob_3p3(double dm, SBN_detector *, int which_dm);


};




#endif
