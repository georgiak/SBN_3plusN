// This macro takes chi2's from selectoins of datasets, applies them to a nice grid and marginalizes  them according to sin22th(mm) or sin22th(em)

#include "fitter.h"

int ntupleProcess(std::string xml){

  bool debug = false;

  ProcessReader rdr;
  if(rdr.Load(xml))
    return 0;

  // Now, sum up the chi2's of everything and throw them in a vector
  float dof;
  std::vector < float > v_chi2, v_mnu, v_um4, v_ue4;

  // Initialize tree variables
  float _chi2, _m_sterile[3], _um_sterile[3], _ue_sterile[3], _dof;

	float chi2_min_total, mnu_min_total, ue_min_total, um_min_total;
	int entry_index_min;

	// Get info from total
	TTree * t_total = (TTree*)rdr.GetFile(0)->Get("Total");
  t_total->SetBranchAddress("chi2",&_chi2);
  t_total->SetBranchAddress("m_sterile",&_m_sterile);
  t_total->SetBranchAddress("um_sterile",&_um_sterile);
  t_total->SetBranchAddress("ue_sterile",&_ue_sterile);
	t_total->SetBranchAddress("dof",&_dof);
	chi2_min_total = 9999;
	for(int j = 0; j < t_total->GetEntries(); j++){
		t_total->GetEntry(j);
		if(_chi2 < chi2_min_total){
			chi2_min_total = _chi2;
			mnu_min_total = _m_sterile[0];
			um_min_total = _um_sterile[0];
			ue_min_total = _ue_sterile[0];
		}
	}
	std::cout << "Best Fit chi2: " << chi2_min_total << ", dof: " << _dof <<  ", mnu: " << mnu_min_total << ", ue: " << ue_min_total << ", um: " << um_min_total << std::endl;

	// Go through all individual experiments
	std::cout << "Time to cycle through the experiments!" << std::endl;
	std::cout << "============================================================" << std::endl;
	float chi2_min, mnu_min, ue_min, um_min;

	for(int i = 0; i < rdr.GetNTrees(); i++){
  	std::cout << rdr.GetName(i) << ":" << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;

    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("m_sterile",&_m_sterile);
    rdr.GetTree(i)->SetBranchAddress("um_sterile",&_um_sterile);
    rdr.GetTree(i)->SetBranchAddress("ue_sterile",&_ue_sterile);
		rdr.GetTree(i)->SetBranchAddress("dof",&_dof);

		// First, find contribution to total chi2 min
		rdr.GetTree(i)->GetEntry(entry_index_min);
		std::cout << "Best fit chi2 contribution: " << _chi2 << std::endl;

		chi2_min = 999999;
    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);
			if(_chi2 < chi2_min){
				chi2_min = _chi2;
				mnu_min = _m_sterile[0];
				um_min = _um_sterile[0];
				ue_min = _ue_sterile[0];
			}
		}
		std::cout << "Dataset chi2 min: " << chi2_min << ", dof: " << _dof <<  ", mnu: " << mnu_min << ", ue: " << ue_min << ", um: " << um_min << std::endl;
		std::cout << "=============================================================" << std::endl;
	}


  return 0;
}

int main(int argc, char* argv[]){

  std::string xml = "";
  int iarg = 0;
  opterr=1;
  int index;

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

  ntupleProcess(xml);
  return 0;
}
