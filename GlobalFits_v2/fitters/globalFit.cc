#include "fitter.h"

int globalFit(std::string xml){

  bool debug = true;

  // Read our XML file
  FitReader rdr;
  if(rdr.Load(xml))
    return 0;
  Oscillator osc = rdr.GetOscillator();

  // Initialize datasets
  std::cout << "Initializing " << rdr.GetNDatasets() << " datasets!" << std::endl;
  int ndf = 0;
  //std::string dataLoc = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/";
  std::string dataLoc = "../../data/";
  for(int i = 0; i < rdr.GetNDatasets(); i++){
    ndf += rdr.GetDataset(i)->Init(dataLoc,osc,debug);
  }
  std::cout << "Datasets Initialized!" << std::endl;

  // Create output File
  std::string outfile = "globFit.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  OutTree * chi2Nt = new OutTree("Total");

  // Create a nu model
  neutrinoModel nuModel, nuModelOld;
  double chi2, chi2Log, chi2LogOld;

  std::cout << "Number of MC models = " << osc.nMCGen  << std::endl;

  std::cout << "Initializing Markov chain parameters..." << std::endl;
  osc.PrintMarkovSettings();
  chi2Log = 0;  chi2LogOld = 0;
  nuModelOld = osc.InitializeMarkovParams();

  std::cout << "Start generating neutrino mass and mixing models" << std::endl;
  int count = 0;
  for(int iMCGen = 1; iMCGen <= osc.nMCGen; iMCGen++){
    std::cout << "Progress: " << iMCGen << " / " << osc.nMCGen << "\r";

    if(count == 0)
      nuModel = osc.InitializeMarkovParams();
    else
      nuModel = osc.NewModel(nuModelOld);

    // Calculate chi2s
    chi2 = 0;
    for(int i = 0; i < rdr.GetNDatasets(); i++){
      chi2 += rdr.GetDataset(i)->Chi2(osc,nuModel,debug);
    }

    chi2Nt->Fill(chi2,ndf,nuModel);
    count ++;
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
	  {0,			no_argument, 		0,  0},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

    switch(iarg){
		  case 'x':
			  xml = optarg;
			  break;
      case 't':
        massStart = atoi(optarg);
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

  globalFit(xml);
  return 0;
}
