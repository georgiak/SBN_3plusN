// fluxprop.cxx by Davio Cianci
// 9/1/2017
//
// This file takes the initial unpropagated flux of either pions or kaons (by using -p or -k in execution) and propagates them through the earth
// Then it takes the propagated flux at icecube and bin averages according to icecube's true neutrino energy resolution
// Then it convolves the resulting flux (in E_nu, CosTheta) with the systematic response array
// The result is an array of proxy E_mu vs CosTheta which can be used to weigh the monte carlo to get the predicted event rate at icecube
//
//
//  Good things to remember. True/Reco Nu Energy Bins: 200, Cos(Theta) bins: 21, Proxy Reco Mu Energy Bins: 10
//
//////////////////////////////////////////////////////////////////////////////////////////

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
  unsigned int numneu = 3;  // number of neutrinos
  bool interactions = true;
  std::string partsrc;
  if (part == 0)
    partsrc = "pion";
  else if (part ==1)
    partsrc = "kaon";


  //
  // First,  load up the initial flux!
  //
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

  //
  // Set up nuSQuiDS
  //
  //Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the
  //units, if they are omitted the input is in eV.
  double Emin=1*units.GeV;
  double Emax=1e6*units.GeV;
  double czmin=-1.;
  double czmax=.24;
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

  // Here is where we'd put in NSI or a 4h neutrino. There are details in the nusquids paper

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
        inistate[ci][ei][0][flv] = (flv == 1) ? vec_nuflux[ei][ci] : 0.0;//set nonzero only to the muon flavor
        inistate[ci][ei][1][flv] = (flv == 1) ? vec_nubarflux[ei][ci] : 0.0;//set nonzero only to the muon flavor
      }
    }
  }

  //Set the initial state in the atmSQuIDS object
  nus_atm.Set_initial_state(inistate,flavor);
  std::cout << "End: setting initial state." << std::endl;

  //Set to true the monitoring prgress bar and the vacuum oscillations
  nus_atm.Set_ProgressBar(false);
  nus_atm.Set_IncludeOscillations(true);

  // Propagate! Everything!
  nus_atm.EvolveState();


  //
  // Get bin averaged, propagated flux at icecube detector
  //
  // We're gonna average 4 evaluations per bin. Energy is true Nu Energy.
  // Number of bins come from resolution of detector, so don't toy too much.
  int nZBins = 21;
  int nEBins = 200;
  int nAvg = 4;           // How many measurements we average together per bin

  // Bin edges
  float binsE[nEBins+1], binsCosZ[nZBins+1];
  binsCosZ[0] = -1.0;
  for(int i = 0; i < nZBins; i++){
    binsCosZ[i+1] = -.96 + i*.06;
  }
  for(int j = 0; j < nEBins+1; j++){
    binsE[j] = pow(10,log10(200) + j*log10(1e6/200)/nEBins);
  }

  std::vector <std::vector < float > > atmFlux_nu, atmFlux_nubar;
  atmFlux_nu.resize(21,std::vector<float>(200));
  atmFlux_nubar.resize(21,std::vector<float>(200));

  float flux_nu, flux_nubar, centerCZ, centerE;
  for(int cz = 0; cz < nZBins; cz++){
    for(int re = 0; re < nEBins; re++){
      flux_nu = 0;  flux_nubar = 0;
      // Average together the contents of the bins...
      for(int i = 0; i < nAvg; i++)
        for(int j = 0; j < nAvg; j++){
          float binCZ = binsCosZ[cz] + i*(binsCosZ[cz+1]-binsCosZ[cz])/nAvg;
          float binE = pow(10,log10(binsE[re]) + j*log10(binsE[re+1]/binsE[re])/nAvg);
          flux_nu += nus_atm.EvalFlavor(1,binCZ, binE*units.GeV, 0);
          flux_nubar += nus_atm.EvalFlavor(1,binCZ, binE*units.GeV, 1);
        }
      // Get energy and cos(zenith) at center of bin
      centerCZ = binsCosZ[cz] + (binsCosZ[cz+1]-binsCosZ[cz])/2;
      centerE  = pow(10,log10(binsE[re]) + log10(binsE[re+1]/binsE[re])/2);
      atmFlux_nu[cz][re] = flux_nu/(nAvg*nAvg);
      atmFlux_nubar[cz][re] = flux_nubar/(nAvg*nAvg);
    }
  }

  //
  //  Not done yet. Now, we invoke the systematic response array!
  //
  std::vector < std::vector < std::vector <float> > > systematicResponseArrayNu, systematicResponseArrayNubar;
  systematicResponseArrayNu.resize(10,std::vector  < std::vector <float> >(21,std::vector <float>(200)));
  systematicResponseArrayNubar.resize(10,std::vector  < std::vector <float> >(21,std::vector <float>(200)));

  f.open("detectorResponseArray.txt");
  // Skip the first 3 lines
  for(int i = 0; i < 3; i++)
    f.ignore(300,'\n');
  float dummy;
  std::string sdummy;
  for(int cz = 0; cz < 21; cz++)  for(int re = 0; re < 200; re++){
    // first, nu
    f >> sdummy;
    f >> dummy;
    f >> dummy;
    for(int pe = 0; pe < 10; pe++)
      f >> systematicResponseArrayNu[pe][cz][re];
    // nubar!
    f >> sdummy;
    f >> dummy;
    f >> dummy;
    for(int pe = 0; pe < 10; pe++)
      f >> systematicResponseArrayNubar[pe][cz][re];
  }
  f.close();

  //
  // Convolve the systematics with the flux
  //
  // This is our final flux. It is in terms of proxy reco mu energy and cos zenith angle
  std::vector < std::vector <float > > flux_conv_nu, flux_conv_nubar;
  flux_conv_nu.resize(21,std::vector<float>(10));
  flux_conv_nubar.resize(21,std::vector<float>(10));


  float eflux_nu, eflux_nubar;
  for(int cz = 0; cz < 21; cz++){
    for(int pe = 0; pe < 10; pe ++){
      eflux_nu = 0;
      eflux_nubar = 0;
      for(int re = 0; re < 200; re ++){
        eflux_nu += systematicResponseArrayNu[pe][cz][re] * atmFlux_nu[cz][re];
        eflux_nubar += systematicResponseArrayNubar[pe][cz][re] *atmFlux_nu[cz][re];
      }
      flux_conv_nu[cz][pe] = eflux_nu;
      flux_conv_nubar[cz][pe] = eflux_nubar;
    }
  }

  //
  //  Add the nu and nubar together and spit out a 21 x 10 array in a nice text file. THEN WE'RE DONE!
  //
  std::ofstream outfile;
  outfile.open ("flux_convolved_"+partsrc+".txt");
  for(int cz = 0; cz < 21; cz++){
    for(int pe = 0; pe < 10; pe++){
      outfile << flux_conv_nu[cz][pe] + flux_conv_nubar[cz][pe];
    }
  }
  outfile.close();

// TODO: throw in the corresponding cosZenith and proxy energy for each element so it's at all usable. now it's just a straight list.


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
