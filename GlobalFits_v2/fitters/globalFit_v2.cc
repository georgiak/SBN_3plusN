#include "MiniBooNE.h"

int main(){

  bool debug = false;

  // Declare datasets
  MiniBooNE mb_nu(false);
  MiniBooNE mb_nubar(true);

  // Initialize datasets
  int ndf = 0;
  std::string dataLoc = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/";

  ndf += mb_nu.Init(dataLoc,debug);
  ndf += mb_nubar.Init(dataLoc,debug);

  // Initialize our oscillator
  Oscillator osc(100.f,.01f,0.f,.5f,.3f,.2f,7.f,1,0,0);

  // Create output File
  std::string outfile = "brute.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  OutTree chi2Nt("Total");

  // Create a nu model
  neutrinoModel nuModel;
  int count = 0;
  float chi2;
  for(int mi = 0; mi < 100; mi++){
    for(int uei = 0; uei < 100; uei++){
      for(int umi = 0; umi < 100; umi++){
        std::cout << "Progress: " << float(count)/(100.*100.) << "\% \r";

		    nuModel.zero();
		    nuModel.Ue[0] = uei/float(100)*(.5);
		    nuModel.Um[0] = umi/float(100)*(.5);
		    nuModel.mNu[0] = pow(10,(mi/float(100)*TMath::Log10(10./.1) + TMath::Log10(.1)));

        // Calculate chi2s
        chi2 = 0;
        chi2 += mb_nu.Chi2(osc,nuModel,debug);
        chi2 += mb_nubar.Chi2(osc,nuModel,debug);

        count ++;
      }
    }
  }

  // Write everything to File
  chi2Nt.Tree()->Write();
  mb_nu.Tree()->Write();
  mb_nubar.Tree()->Write();
  f->Close();

  return 0;
}
