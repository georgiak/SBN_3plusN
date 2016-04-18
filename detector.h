#ifndef DETECTOR_H_
#define DETECTOR_H_

#include "TRandom.h"
#include <vector>
#include <string>

/*************************************************************
 *************************************************************
 *	TODO:
 *	    (4) overload detectors so I can just pass identifiers DONE
 ************************************************************
 ************************************************************/

#define DET_SBND 0
#define DET_UBOONE 1
#define DET_ICARUS 2

	double smear_energy(double En, double Percen, TRandom * rangen);
	double muon_track_length(double El);
	double pion_track_length(double El);
	double photon_conversion_length(double ep, TRandom * r);
	double bethe(double beta);

	double CSDA_integrand(double Emu);
	double pion_containment(double posX, double posY, double posZ, TRandom * r);

	int get_endpoint(double *vertex,double track_L,double * pl,double *  endpoint);



class SBN_detector {
	double dh, dw, dl;


	public:
	char const * name;
	bool mumode;
	double height, width, length;
	double f_height, f_width, f_length;
        double volume;
	double f_volume;	
	double baseline;
	char const * fname; // location of root ntuple containing data file		
	char const * foscname;	//location of full oscillated 
	double potmodifier;
	int identifier;

	//Constructors, if blank corresponds to background only. 
	SBN_detector (double h, double w, double l, double fh, double fw, double fl,double base);
	SBN_detector (int identifier, bool ismue = false);

	double osc_length(TRandom * rangen);

	bool is_active(double * pos);
	bool is_fiducial(double * pos);
	void random_pos(TRandom * rangen, double * vec);

	double track_length_escape(double * inside, double * outside);
	bool is_fully_contained(double *vertex,double * endpoint);
	


};



#endif
