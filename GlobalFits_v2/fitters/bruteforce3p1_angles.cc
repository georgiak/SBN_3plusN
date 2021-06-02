 #include "fitter.h"

bool debug;;

int bruteforce(std::string xml){

  // Read our XML file
  FitReader rdr;
  if(rdr.Load(xml))
    return 0;
  Oscillator osc = rdr.GetOscillator();

  // Initialize datasets
  std::cout << "Initializing " << rdr.GetNDatasets() << " datasets!" << std::endl;
  int ndf = 0;
  std::string dataLoc = "/home/dcianci/Physics/GlobFitDocumentation/SBN_3plusN/GlobalFits_v2/";
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
  neutrinoModel nuModel;

  int count = 0;
  int grdpts = osc.GridSize();
  float chi2;
  std::cout << "Beginning chi2 loop" << std::endl;

  double theta_lowbound(.01), theta_hibound(3.1415926/4.);
  double mnu_lowbound(.1), mnu_hibound(10.00);

  for(int mi = 0; mi < grdpts; mi++){
    for(int t41i = 0; t41i < grdpts; t41i++){
      for(int t42i = 0; t42i < grdpts; t42i++){
        if(!debug)
          std::cout << "Progress: " << float(count)/(pow(grdpts,3)/(100.f)) << "\% \r";

  	    double mnu = IndexToValue(mi,mnu_lowbound,mnu_hibound,grdpts);
        double theta14 = IndexToValue(t41i,theta_lowbound,theta_hibound,grdpts);
        double theta24 = IndexToValue(t42i,theta_lowbound,theta_hibound,grdpts);
        double theta34 = 0;
        nuModel.Init(mnu,theta14,theta24,theta34);

        chi2 = 0;
        // Calculate chi2s
        for(int i = 0; i < rdr.GetNDatasets(); i++){
          chi2 += rdr.GetDataset(i)->Chi2(osc,nuModel,debug);
        }

        if(debug){
          std::cout << chi2 << std::endl;
        }

        chi2Nt->Fill(chi2,ndf,nuModel);
        count ++;
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

  bruteforce(xml);
  return 0;
}
