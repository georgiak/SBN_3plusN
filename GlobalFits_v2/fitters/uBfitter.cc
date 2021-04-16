#include "fitter.h"
#include "MicroBooNE_dis.h"

bool debug;

int bruteforce(){

  debug = true;
  const int ngrdpts = 500;
  MicroBooNE_dis* ub = new MicroBooNE_dis;

  Oscillator osc;
  osc.gridpts = ngrdpts;

  std::cout << "Initializing"  << std::endl;
  std::string dataLoc = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/";
  ub->Init(dataLoc,osc,debug);
  std::cout << "Dataset Initialized!" << std::endl;

  // Create output File
  std::string outfile = "ubfit.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  OutTree * chi2Nt = new OutTree("Total");

  // Create a nu model
  neutrinoModel nuModel;

  int count = 0;
  float chi2;
  std::cout << "Beginning chi2 loop" << std::endl;

  double theta_lowbound(.01), theta_hibound(3.1415926/4.);
  double mnu_lowbound(.1), mnu_hibound(1000.00);

  for(int mi = 0; mi < ngrdpts; mi++){
      for(int t42i = 0; t42i < ngrdpts; t42i++){
        if(!debug)
          std::cout << "Progress: " << float(count)/(pow(ngrdpts,2)/(100.f)) << "\% \r";
        mi = 481;
  	    double mnu = pow(10,(mi/float(ngrdpts)*TMath::Log10(mnu_hibound/mnu_lowbound) + TMath::Log10(mnu_lowbound)));
        double theta14 = 0;
        double theta34 = 0;
        double sin22th = 0.8;//pow(10,(t42i/float(ngrdpts)*TMath::Log10(1.0/0.01) + TMath::Log10(0.01)));
        double theta24 = asin(sqrt((1.0-sqrt(1.0-sin22th))/2.0));
        nuModel.Init(mnu,theta14,theta24,theta34);

        chi2 = 0;
        // Calculate chi2s
        chi2 = ub->Chi2(osc,nuModel,debug);

        if(debug){
          //std::cout << "Total chi2 for model: " << mnu  << " " <<  theta14 << " " << theta24 << " " <<  theta34 << std::endl;
          std::cout << chi2 << std::endl;
        }

        chi2Nt->Fill(chi2,0,nuModel);
        count ++;
        return 1;
      }
  }

  // Write everything to File
  std::cout << "Writing to file." << std::endl;
  chi2Nt->Write();
  ub->Write();
  f->Close();

  return 0;
}

int main(){

  bruteforce();
  return 0;
}
