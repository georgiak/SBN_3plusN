#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "nuSQuIDS/nuSQuIDS.h"

using namespace nusquids;

int fluxprop(int part){

  squids::Const units;
  unsigned int numneu = 3;
  bool interactions = true;
  std::string partsrc;
  if (part == 0)
    partsrc = "pion";
  else if (part ==1)
    partsrc = "kaon";


  std::vector < std::vector < float > > vec_energy, vec_zenith, vec_nuflux, vec_nubarflux;
  vec_energy.resize(150,std::vector<float>(40));
  vec_zenith.resize(150,std::vector<float>(40));
  vec_nuflux.resize(150,std::vector<float>(40));
  vec_nubarflux.resize(150,std::vector<float>(40));

  std::ifstream f;
  f.open("IC86SterileNeutrinoDataRelease/atmospheric_flux/initial/initial_"+partsrc+"_atmopheric_HondaGaisser.dat");
  for(int cosz = 0; cosz < 40; cosz++){
    for(int recoe = 0; recoe < 150; recoe++){
      f >> vec_zenith[recoe][cosz];
      f >> vec_energy[recoe][cosz];
      f >> vec_nuflux[recoe][cosz];
      f >> vec_nubarflux[recoe][cosz];
    }
  }
  f.close();

  //Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the
  //units, if they are omitted the input is in eV.
  double Emin=1*units.GeV;
  double Emax=1e6*units.GeV;
  double czmin=-1.;
  double czmax=.2;
  //Declaration of the atmospheric object
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,40),logspace(Emin,Emax,150),numneu,both,interactions);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;

  std::cout << "Begin: setting mixing angles." << std::endl;
  // set mixing angles, mass differences and cp phases
  nus_atm.Set_MixingAngle(0,1,0.563942);
  nus_atm.Set_MixingAngle(0,2,0.154085);
  nus_atm.Set_MixingAngle(1,2,0.785398);

  nus_atm.Set_SquareMassDifference(1,7.65e-05);
  nus_atm.Set_SquareMassDifference(2,0.00247);

  nus_atm.Set_CPPhase(0,2,0);

  //Setup integration precision
  nus_atm.Set_rel_error(1.0e-6);
  nus_atm.Set_abs_error(1.0e-6);
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Array that contains the values of the energies and cosine of the zenith, is the same length for every zenith
  auto e_range = nus_atm.GetERange();
  auto cz_range = nus_atm.GetCosthRange();

  std::cout << "Begin: setting initial state." << std::endl;

  //Construct the initial state, we set a flat spectra in zenith and log-energy
  marray<double,4> inistate{nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);
  for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
      for (int flv = 0; flv < numneu; flv++){
        inistate[ci][ei][0][flv] = (flv == 1) ? vec_nuflux[ei][ci] : 0.0;//set 1 only to the muon flavor
        inistate[ci][ei][1][flv] = (flv == 1) ? vec_nubarflux[ei][ci] : 0.0;//set 1 only to the muon flavor
      }
    }
  }

  //Set the initial state in the atmSQuIDS object
  nus_atm.Set_initial_state(inistate,flavor);
  std::cout << "End: setting initial state." << std::endl;

  //Set to true the monitoring prgress bar and the vacuum oscillations
  nus_atm.Set_ProgressBar(false);
  nus_atm.Set_IncludeOscillations(true);

/*
  // Here's our sterile boy
  for(int im = 0; im < 100; im++) for(int is = 0; is < 100; is++){

    float sin22th = pow(10,(float(is)/100.*log10(1./.01) + log10(.01)));
    float dm2 = pow(10,(float(im)/100.*log10(100./.01) + log10(.01)));

    std::cout << "sin22th: " << sin22th << " dm2: " << dm2  <<  std::endl;
    nus_atm.Set_SquareMassDifference(3,dm2);
    nus_atm.Set_MixingAngle(1,3,sin22th);
*/

    //Set the initial state in the atmSQuIDS object
    nus_atm.Set_initial_state(inistate,flavor);
    nus_atm.EvolveState();

    //This file will contain the final flux, since initially we set it to 1 for the muon,
    //this can be read as the muon ration F_final/F_initial in cos(zenith) and energy.
    //std::ofstream file("libflux/flux_"+partsrc+"_dm2_"+std::to_string(im)+"_sin22th_"+std::to_string(is)+".txt");
    std::ofstream file("libflux/flux_"+partsrc+"_nominal.txt");

    //Set the resolution and the ranges for the ouput, remember that an interpolation in energy is used in
    //in the interaction picture, the vacuum oscillations are solve analytically with arbitrary Energy precision.
    //For the zenith a linear interpolation is used.
    int Nen=700;
    int Ncz=100;
    double lEmin=log10(Emin);
    double lEmax=log10(Emax);;

    //Writing to the file!
    file << "# cos(zenith) E flux_nu flux nubar . . . ." << std::endl;
    for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
      for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
        double E=pow(10.0,lE);
        file << " " << cz << " " << E/units.GeV;
        file << " " <<  nus_atm.EvalFlavor(1,cz, E, 0);
        file << " " <<  nus_atm.EvalFlavor(1,cz, E, 1);
        file << std::endl;
      }
    }
    file.close();

  return 0;
}


int main(int argc, char* argv[])
{
  int part, index;
  part = -1;
  int iarg = 0;
  opterr=1;
  const struct option longopts[] = {
    {"pion",	 		no_argument, 		0, 'p'},
    {"kaon",	 	no_argument, 		0, 'k'},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "pk", longopts, &index);

    switch(iarg)
    {
      case 'p':
        part = 0;
        break;
      case 'k':
        part = 1;
        break;
    }
  }
  if(part == -1){
    std::cout << "Pion (-p)  or kaon (-k), my man?" << std::endl;
    return 0;
  }
  else std::cout << part << std::endl;
  fluxprop(part);
  return 0;
}
