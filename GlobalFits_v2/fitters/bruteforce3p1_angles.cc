 #include "fitter.h"

bool debug;;

int bruteforce(std::string xml, int massStart = -1){

  // Read our XML file
  FitReader rdr;
  if(rdr.Load(xml))
    return 0;
  Oscillator osc = rdr.GetOscillator();

  // Initialize datasets
  std::cout << "Initializing " << rdr.GetNDatasets() << " datasets!" << std::endl;
  int ndf = 0;
  //std::string dataLoc = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/";
  std::string dataLoc = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/";
  for(int i = 0; i < rdr.GetNDatasets(); i++){
    ndf += rdr.GetDataset(i)->Init(dataLoc,osc,debug);
  }
  std::cout << "Datasets Initialized!" << std::endl;

  // Create output File
  std::string outfile = "brute.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  OutTree * chi2Nt = new OutTree("Total");

  // Create a nu model
  //neutrinoModel nuModel;

  int count = 0;
  int grdpts = osc.GridSize();
  float chi2;
  std::cout << "Beginning chi2 loop" << std::endl;
  int mStart, mEnd;
  if(massStart < 0){
    mStart = 0; mEnd = grdpts;
  }

  double thetamin(1e-4), thetamax(3.14159/2);
  //double thetamin(3e-2), thetamax(3e-1);

  for(int mi = mStart; mi < mEnd; mi++){
    for(int t41i = 0; t41i < grdpts; t41i++){
      for(int t42i = 0; t42i < grdpts; t42i++){

        std::cout << "Progress: " << float(count)/(pow(grdpts,3)/(100.f)) << "\% \r";

        neutrinoModel nuModel;
		    double mnu = pow(10,(mi/float(grdpts)*TMath::Log10(10./.1) + TMath::Log10(.1)));
        // new min for mnu
        double theta14 = pow(10,(t41i/float(grdpts)*TMath::Log10(thetamax/thetamin) + TMath::Log10(thetamin)));
        double theta24 = pow(10,(t42i/float(grdpts)*TMath::Log10(thetamax/thetamin) + TMath::Log10(thetamin)));
        double theta34 = 0;
        nuModel.Init(mnu,theta14,theta24,theta34);
        //nuModel.Init(1.0, 0, 0,0.0);

        //if(osc.RejectModel(nuModel))
        //  continue;

        // Calculate chi2s
        chi2 = 0;
        for(int i = 0; i < rdr.GetNDatasets(); i++){
          chi2 += rdr.GetDataset(i)->Chi2(osc,nuModel,debug);
        }

        chi2Nt->Fill(chi2,ndf,nuModel);
        count ++;
        //return 1;
      }
    }
  }

  // Write everything to File
  std::cout << "Writing to file." << std::endl;
  chi2Nt->Write();

  for(int i = 0; i < rdr.GetNDatasets(); i++){
    rdr.GetDataset(i)->Write();
  }

  f->Close();

  return 0;
}

int main(int argc, char* argv[]){

  std::string xml = "";
  int iarg = 0;
  opterr=1;
  int index;
  int massStart = -1;

  const struct option longopts[] = {
    {"xml", 		required_argument, 	0, 'x'},
    {"debug", optional_argument, 0, 'd'},
	  {0,			no_argument, 		0,  0},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "x:t:d", longopts, &index);

    switch(iarg){
		  case 'x':
			  xml = optarg;
			  break;
      case 't':
        massStart = atoi(optarg);
        break;
      case 'd':
        debug = true;
        std::cout << "DEBUG MODE" << std::endl;
        break;
      case '?':
		  case 'h':
			  std::cout<<"I need an input, friend."<<std::endl;
			  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
			  return 0;
	  }
  }
  if(xml == ""){
    std::cout << "Gimme an XML input or I won't start, I swear to god." << std::endl;
    return 0;
  }

  bruteforce(xml,massStart);
  return 0;
}
