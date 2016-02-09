#ifndef MODEL_H_
#define MODEL_H_

#include <cmath>
#include <vector>
#include "TH1F.h"
#include <iostream>


/*********************************************
* Data Struct: Same as Davio's output.
* -------------------------------------------
* ToDo: 
* Should be packed in a nTuple to me later. 
* 
* *******************************************/


struct neutrinoModel{
	double mNu[3], Ue[3], Um[3], phi[3];
	double dm41Sq, dm51Sq, dm61Sq, dm54Sq, dm64Sq, dm65Sq;
	void zero(){
		for(int i = 0; i < 3; i ++){
			mNu[i] = 0; Ue[i] = 0;
			Um[i] = 0;  phi[i] = 0;
		}
	}
	void difference(){
		dm41Sq = pow(mNu[0],2);
		dm51Sq = pow(mNu[1],2);
		dm61Sq = pow(mNu[2],2);
		dm54Sq = dm51Sq - dm41Sq;
		dm64Sq = dm61Sq - dm41Sq;
		dm65Sq = dm61Sq - dm51Sq;
	}
};

neutrinoModel initializeMarkovParams();
neutrinoModel newModel(neutrinoModel modelOld);
bool rejectModel(neutrinoModel model);

/*********************************************
* Main Class  SBN_spectrum
* -------------------------------------------
* ToDo:
* 
* 
* *******************************************/

class SBN_spectrum {
	struct neutrinoModel nullModel;
	static const int N_m_bins = 19;
	static const int N_e_bins = 11;
	static constexpr double mu_bins[N_m_bins+1] =  {0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.25,1.5,2,2.5};
	static constexpr double e_bins[N_e_bins+1]= {0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.3,1.5,1.75,2,3};

	public:
	std::vector<double > muboone_e 	,sbnd_e, icarus_e , sbnd_m ,  muboone_m , icarus_m ;
	TH1D * THmuboone_e,* THmuboone_m, * THsbnd_e, * THsbnd_m,* THicarus_e,* THicarus_m;


	struct neutrinoModel workingModel;	
	
	//Constructors, if blank corresponds to background only. 
	SBN_spectrum (struct neutrinoModel nuModel);
	SBN_spectrum ();

	std::vector<double > get_sixvector();
	void update_model(struct neutrinoModel nuModel);
	
	//ToDo
	void oscillate();
	void fill_hists();
        void fill_vectors();	
	void THprint();

};

#endif
