// This macro takes chi2's from selectoins of datasets, applies them to a nice grid and marginalizes  them according to sin22th(mm) or sin22th(em)

#include "fitter.h"
#include "TH3D.h"

int ntupleProcess(){

  bool debug = false;

  std::array<float,240000> chi2vec, dm2vec, t14vec, t24vec, t34vec;

  std::string s_dummy,line;
  double d_dummy;
  double loglmin(9999999), odm2min,ot14min,ot24min;
  // okay. load up txt  file.
  ifstream file;
  file.open("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/logl_final_hopefully_nov30.txt");
  for(int i = 0; i < 240000; i++){
    file >> s_dummy;
    file >> dm2vec[i];
    file >> t14vec[i];
    file >> t24vec[i];
    file >> t34vec[i];
    file >> d_dummy;
    file >> d_dummy;
    file >> chi2vec[i];
    if(chi2vec[i] < loglmin){
      loglmin = chi2vec[i];
      odm2min = dm2vec[i];
      ot14min = t14vec[i];
      ot24min = t24vec[i];
    }
  }
  file.close();

  std::cout << "LOGL MIN: " << loglmin << std::endl;
  std::cout << "(" << odm2min << ", " << ot14min << ", " << ot24min << ")" << std::endl;

  return 1;

  if(debug == true){

    float mstep = TMath::Log10(100.f/.01f)/25;
    float sin22step = TMath::Log10(1.f/float(1e-5))/25;
    std::vector < std::vector < float > > chi2grid;
    chi2grid.assign(25, std::vector < float > (25, 0.));

    for(int i = 0; i < chi2vec.size(); i++){
      double mysin22th = 4 * pow(cos(t14vec[i]),2) * pow(sin(t24vec[i]),2) * (1 - pow(cos(t14vec[i]),2) * pow(sin(t24vec[i]),2));
      //if(i%20==0) std::cout << mysin22th << std::endl;
      int is = (int)TMath::Nint(TMath::Log10(mysin22th/1e-5)/sin22step);
      int im = (int)TMath::Nint(TMath::Log10(dm2vec[i]/.01)/mstep);

      if(is < 25 && im < 25 && is > -1){
        if(chi2grid[is][im] == 0)
          chi2grid[is][im] = chi2vec[i];
        else
          chi2grid[is][im] = TMath::Min(chi2grid[is][im],chi2vec[i]);
      }
      else
        continue;
    }

    // Print out chiogram
    for(int i = 0; i < 25; i++){
      for(int j = 0; j< 25; j++){
        double _sin22th = pow(10,(i/float(25) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
        double _dm2 = pow(10,(j/float(25) * TMath::Log10(100./.01) + TMath::Log10(.01)));
        std::cout << "B: " << _sin22th << " " << _dm2 << " " << chi2grid[i][j] << std::endl;
      }
    }
  }

  // Ok. very rad. i... think we're done here?
  // Create output File
  std::string outfile = "IC_chi2grids.root";
	//std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  // Cool. Now, let's make an array of 3d histograms: one for each theta34. These histograms will be in log10([dim]) so linear interpolation will work.
  std::vector<TH3D*> vecLogTHisto;

  // define bin edges for the histograms.
  const double logm41bins[41] = {-1.  , -0.95, -0.9 , -0.85, -0.8 , -0.75, -0.7 , -0.65, -0.6 ,
                                -0.55, -0.5 , -0.45, -0.4 , -0.35, -0.3 , -0.25, -0.2 , -0.15,
                                -0.1 , -0.05,  0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,
                                0.35,  0.4 ,  0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,
                                0.8 ,  0.85,  0.9 ,  0.95,  1.  };




  //  for thetas 14 and 24, you'll notice that the beginning and end don't match exactly with the theta bounds.
  //  That's on purpose. I added one buffer bin on either side for easier interpolation.
  //  The second entry from either end are the actual theta bounds
  //  I... forgot to do this for m41, but it's less important there since we're not as concerned with super high or super low mass
  //  I'd redo it but it takes literally five days to run.
  //  And, of  course, there's no interpolation with theta34

  const double logtheta24bins[41] = {-2.04987079, -2.        , -1.95012921, -1.90025843, -1.85038764,
                                    -1.80051686, -1.75064607, -1.70077528, -1.6509045 , -1.60103371,
                                    -1.55116292, -1.50129214, -1.45142135, -1.40155057, -1.35167978,
                                    -1.30180899, -1.25193821, -1.20206742, -1.15219664, -1.10232585,
                                    -1.05245506, -1.00258428, -0.95271349, -0.9028427 , -0.85297192,
                                    -0.80310113, -0.75323035, -0.70335956, -0.65348877, -0.60361799,
                                    -0.5537472 , -0.50387642, -0.45400563, -0.40413484, -0.35426406,
                                    -0.30439327, -0.25452248, -0.2046517 , -0.15478091, -0.10491013,
                                    -0.05503934};

  const double logtheta14bins[16] = {-2.14577614, -2.        , -1.85422386, -1.70844771, -1.56267157,
                                    -1.41689542, -1.27111928, -1.12534314, -0.97956699, -0.83379085,
                                    -0.6880147 , -0.54223856, -0.39646241, -0.25068627, -0.10491013,
                                    0.04086602};

  double m41min(pow(10,logm41bins[0])), m41max(pow(10,logm41bins[40]));
  double t14min(pow(10,logtheta14bins[0])), t14max(pow(10,logtheta14bins[15]));
  double t24min(pow(10,logtheta24bins[0])), t24max(pow(10,logtheta24bins[40]));
  double t34min(.01), t34max(3.1415926/4.0);

  std::cout << "t14maxmin: " << t14max << " " << t14min << std::endl;
  std::cout << "t24maxmin: " << t24max << " " << t24min << std::endl;
  std::cout << "t34maxmin: " << t34max << " " << t34min << std::endl;


  std::array<std::array<std::array<std::array<float,40>,15>,40>,10> db_vecLogTHisto; // for debug

  for(int i = 0; i < 10; i++){
    vecLogTHisto.push_back(new TH3D(("loghisto_"+std::to_string(i)).c_str(),("loghisto_"+std::to_string(i)).c_str(),40,logm41bins,15,logtheta14bins,40,logtheta24bins));
  }

  // Okay. Now, let's fill all those histograms.
  // Loop through our array
  int im41, ith14, ith24, ith34;
  for(int i = 0; i < 240000; i++){
    im41 = floor(40.0*log10(sqrt(dm2vec[i])/m41min)/log10(m41max/m41min)+1.0);
    ith14 = floor(15*log10(t14vec[i]/t14min)/log10(t14max/t14min)+1.0);
    ith24 = floor(40*log10(t24vec[i]/t24min)/log10(t24max/t24min)+1.0);
    ith34 = floor(10*log10(t34vec[i]/t34min)/log10(t34max/t34min));

    if(im41 == 1 and ith14 == 1 and ith24 == 1){
      std::cout << "F: " << ith34 << " "  << t34vec[i] << std::endl;
    }

    // importantly, we need to convert log likelihood to chi2.
    vecLogTHisto[ith34]->SetBinContent(im41,ith14,ith24,2*(chi2vec[i]-loglmin));
    db_vecLogTHisto[ith34][im41-1][ith14-1][ith24-1] = 2*(chi2vec[i]-loglmin);
  }

  if(debug==true){

    /*
    float mstep = TMath::Log10(100.f/.01f)/25;
    float sin22step = TMath::Log10(1.f/float(1e-5))/25;
    std::vector < std::vector < float > > chi2gridint;
    chi2gridint.assign(25, std::vector < float > (25, 0.));

    for(int i = 0; i < chi2vec.size(); i++){
      double mysin22th = 4 * pow(cos(t14vec[i]),2) * pow(sin(t24vec[i]),2) * (1 - pow(cos(t14vec[i]),2) * pow(sin(t24vec[i]),2));
      //if(i%20==0) std::cout << mysin22th << std::endl;
      int is = (int)TMath::Nint(TMath::Log10(mysin22th/1e-5)/sin22step);
      int im = (int)TMath::Nint(TMath::Log10(dm2vec[i]/.01)/mstep);
      int ith34 = floor(10*log10(t34vec[i]/t34min)/log10(t34max/t34min));

      if(is < 25 && im < 25 && is > -1){
        double chi2int = vecLogTHisto[ith34]->Interpolate(log10(dm2vec[i]),log10(t14vec[i]),log10(t24vec[i]));//min(vecLogICHisto[i]->Interpolate(log10(ops[0]),log10(ops[1]),log10(ops[2])),chi2);)
        if(chi2gridint[is][im] == 0)
          chi2gridint[is][im] = chi2int;
        else
          chi2gridint[is][im] = TMath::Min(chi2gridint[is][im],chi2int);
      }
      else
        continue;
    }

    // Print out chiogram
    for(int i = 0; i < 25; i++){
      for(int j = 0; j< 25; j++){
        double _sin22th = pow(10,(i/float(25) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
        double _dm2 = pow(10,(j/float(25) * TMath::Log10(100./.01) + TMath::Log10(.01)));
        std::cout << "R: " << _sin22th << " " << _dm2 << " " << chi2gridint[i][j] << std::endl;
      }
    }
    */

    // Print out debug histograms
    for(int it34 = 0; it34<10; it34++){
      std::cout << "theta34: " << it34 << std::endl;
      for(int idm2 = 0; idm2<40; idm2++){
        std::cout << "dm2: " << idm2 << std::endl;
        for(int it14 = 0; it14<15; it14++){
          for(int it24 = 0; it24<40; it24++){
            std::cout << db_vecLogTHisto[it34][idm2][it14][it24] << " ";
          }
          std::cout << endl;
        }
      }
    }
  }


  for(int i = 0; i < 10; i++){
    vecLogTHisto[i]->Write();
  }
  f->Close();

  return 0;
}

int main(){
  ntupleProcess();
  return 0;
}
