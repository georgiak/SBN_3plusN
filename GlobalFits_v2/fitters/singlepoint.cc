#include "fitter.h"

int bruteforce(std::string xml, int massStart = -1){

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

  OutTree * chi2Nt = new OutTree("Total");

  // Define point
  neutrinoModel nuModel;
  nuModel.zero();
  nuModel.Ue[0] = sqrt(.96/4.f);
  nuModel.Um[0] = 1;
  nuModel.mNu[0] = sqrt(.04);

  // Calculate chi2s
  float chi2 = 0;
  for(int i = 0; i < rdr.GetNDatasets(); i++){
    chi2 += rdr.GetDataset(i)->Chi2(osc,nuModel,debug);
  }
  chi2Nt->Fill(chi2,ndf,nuModel);

  std::cout << "CHI2: " << chi2 << std::endl;

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

  bruteforce(xml,massStart);
  return 0;
}
