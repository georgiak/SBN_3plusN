// This macro takes chi2's from selectoins of datasets, applies them to a nice grid and marginalizes  them according to sin22th(mm) or sin22th(em)

#include "fitter.h"

bool stringCompare(const string & l, const string & r){
	// We only want to compare the points, not the chisq
	return (l.substr(0,6) == r.substr(0,6));
}

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
  float chi2_min = 99999;

  std::cout << "Combining the chi2s" << std::endl;
  for(int i = 0; i < rdr.GetNTrees(); i++){
    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("m_sterile",&_m_sterile);
    rdr.GetTree(i)->SetBranchAddress("um_sterile",&_um_sterile);
    rdr.GetTree(i)->SetBranchAddress("ue_sterile",&_ue_sterile);
    rdr.GetTree(i)->SetBranchAddress("dof",&_dof);
    if(i == 0){
      v_chi2.resize(rdr.GetTree(i)->GetEntries());
      v_mnu.resize(rdr.GetTree(i)->GetEntries());
      v_um4.resize(rdr.GetTree(i)->GetEntries());
      v_ue4.resize(rdr.GetTree(i)->GetEntries());
    }

    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);

      if(i == 0){
        if(j == 0)  dof = _dof;
        v_chi2[j] = _chi2;
        v_mnu[j] = _m_sterile[0];
        v_um4[j] = _um_sterile[0];
        v_ue4[j] = _ue_sterile[0];
      }
      else{
        if(j == 0)  dof += _dof;
          v_chi2[j] += _chi2;
        if(v_mnu[j] != _m_sterile[0] || v_um4[j] != _um_sterile[0]  || v_ue4[j] != _ue_sterile[0]){
          std::cout << "BOOM " << j << " of " << rdr.GetTree(i)->GetEntries() << std::endl;
          std::cout << v_mnu[j] << " " << _m_sterile[0] << " --- " <<  v_um4[j] << " " << _um_sterile[0] << " --- " << v_ue4[j] << " " << _ue_sterile[0] << std::endl;
          std::cout << "Grids are NOT aligned and this code is not prepared to sort them out." << std::endl;
          return 0;
        }
      }
      chi2_min = min(chi2_min, v_chi2[j]);
    }
  }
  std::cout << "Chi2 min: " << chi2_min << std::endl;

  // Create output File
  std::string outfile = rdr.tag + "_proc.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  // Create our trees
  TNtuple * chi2_99 = new TNtuple("chi2_99","chi2_99","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
  TNtuple * chi2_90 = new TNtuple("chi2_90","chi2_90","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
  chi2_min = 99999;

  // Now, we're going to convert to a string of coordinates.
  std::cout << "Convert entries to strings of coords" << std::endl;

  float mstep = TMath::Log10(10.f/.1f)/float(rdr.gridpts_dm2);
  float ustep = (.5)/float(rdr.gridpts_sin22th);

  std::vector<std::string> coords;
  coords.resize(v_chi2.size());

  int _m4, _ue4, _um4;
  std::string _m4s, _ue4s, _um4s, _chi2s;

  for(int i = 0; i < v_chi2.size(); i++){

    // Now, the chi2!
    stringstream stream;
    stream << fixed << setprecision(2) << v_chi2[i];
    _chi2s = stream.str();
    _m4 = TMath::Max(int(floor(TMath::Log10(v_mnu[i]/.1)/mstep)),0);
    if(_m4 < 10) _m4s = "0" + to_string(_m4);	else _m4s = to_string(_m4);
    _ue4 = TMath::Max(int(floor(v_ue4[i]/ustep)),0);
    if(_ue4 < 10) _ue4s = "0" + to_string(_ue4);	else _ue4s = to_string(_ue4);
    _um4 = TMath::Max(int(floor(v_um4[i]/ustep)),0);
    if(_um4 < 10) _um4s = "0" + to_string(_um4);	else _um4s = to_string(_um4);
    coords[i] = _m4s + _ue4s + _um4s + _chi2s;
    std::cout << float(i)/float(v_chi2.size()) << "\r";
  }

  std::cout << "\nCoords all full up.\nWe've got " <<  coords.size() << " of 'em" << std::endl;
  sort(coords.begin(),coords.end());
  coords.erase(unique(coords.begin(), coords.end(),stringCompare), coords.end());

  std::cout << "Holy shit, that worked.\nNow, fill the new ntuple" << std::endl;
  for(int i = 0; i < coords.size(); i++){
    // Gotta convert those boys from string coordinates into actual positions in phase space
    _m4 = stoi(coords[i].substr(0,2));
    float m4p = pow(10,(_m4/float(rdr.gridpts_dm2)*TMath::Log10(10.f/.1f) + TMath::Log10(.1f)));

    _ue4 = stoi(coords[i].substr(2,2));
    float ue4p = _ue4/float(rdr.gridpts_sin22th)*(.5);

    _um4 = stoi(coords[i].substr(4,2));
    float um4p = _um4/float(rdr.gridpts_sin22th)*(.5);

    float chi2p = stof(coords[i].substr(6,6));

    if(chi2p < chi2_min + 6.25)
      chi2_90->Fill(chi2p,m4p,ue4p,um4p,0,0,0,0,0,0,0,0,0);
    if(chi2p < chi2_min + 11.34)
      chi2_99->Fill(chi2p,m4p,ue4p,um4p,0,0,0,0,0,0,0,0,0);
  }

  std::cout << "For 99%%, we have " << chi2_99->GetEntries() << " entries." <<  std::endl;
  std::cout << "For 90%%, we have " << chi2_90->GetEntries() << " entries." <<  std::endl;

  chi2_90->Write();
  chi2_99->Write();

  f->Close();

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
