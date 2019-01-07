// This macro takes chi2's from selectoins of datasets, applies them to a nice grid and marginalizes  them according to sin22th(mm) or sin22th(em)

#include "fitter.h"

int ntupleProcess(std::string xml){

  bool debug = false;

  ProcessReader rdr;
  if(rdr.Load(xml))
    return 0;

  // Now, sum up the chi2's of everything and throw them in a vector
  float dof;
  std::vector < float > v_chi2, v_mnu4, v_mnu5, v_mnu6, v_um4, v_ue4;

  // Initialize tree variables
  float _chi2, _m_sterile[3], _um_sterile[3], _ue_sterile[3], _dof;
	float _chi2_min(99999.f), _m4_min, _m5_min, _m6_min, _ue4_min, _ue5_min, _ue6_min, _um4_min, _um5_min, _um6_min;

  std::cout << "Combining the chi2s" << std::endl;
  for(int i = 0; i < rdr.GetNTrees(); i++){
		std::cout << "A" << std::endl;
    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("m_sterile",&_m_sterile);
    rdr.GetTree(i)->SetBranchAddress("um_sterile",&_um_sterile);
    rdr.GetTree(i)->SetBranchAddress("ue_sterile",&_ue_sterile);
    rdr.GetTree(i)->SetBranchAddress("dof",&_dof);
    if(i == 0){
      v_chi2.resize(rdr.GetTree(i)->GetEntries());
      v_mnu4.resize(rdr.GetTree(i)->GetEntries());
      v_mnu5.resize(rdr.GetTree(i)->GetEntries());
      v_mnu6.resize(rdr.GetTree(i)->GetEntries());
      v_um4.resize(rdr.GetTree(i)->GetEntries());
      v_ue4.resize(rdr.GetTree(i)->GetEntries());
    }
    rdr.GetTree(i)->GetEntry(90);
		std::cout << "B" << std::endl;
    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);

      if(i == 0){
        if(j == 0)  dof = _dof;
        v_chi2[j] = _chi2;
        v_mnu4[j] = _m_sterile[0];
        v_mnu5[j] = _m_sterile[1];
        v_mnu6[j] = _m_sterile[2];
        v_um4[j] = _um_sterile[0];
        v_ue4[j] = _ue_sterile[0];
				if(v_mnu5[j] < .1 || v_mnu6[j] < .1)
					std::cout << "FUCK: " << v_mnu5[j] << " "  << v_mnu6[j] << std::endl;
      }
      else{
          if(j == 0)  dof += _dof;
          v_chi2[j] += _chi2;
          if(v_mnu4[j] != _m_sterile[0] || v_um4[j] != _um_sterile[0]  || v_ue4[j] != _ue_sterile[0]){
            std::cout << "BOOM " << j << " of " << rdr.GetTree(i)->GetEntries() << std::endl;
            std::cout << v_mnu4[j] << " " << _m_sterile[0] << " --- " <<  v_um4[j] << " " << _um_sterile[0] << " --- " << v_ue4[j] << " " << _ue_sterile[0] << std::endl;
            std::cout << "Grids are NOT aligned and this code is not prepared to sort them out." << std::endl;
            return 0;
          }
      }

			if(_chi2 < _chi2_min){
				_chi2_min = _chi2;
				_m4_min = _m_sterile[0];
				_m5_min = _m_sterile[1];
				_m6_min = _m_sterile[2];
				_ue4_min = _ue_sterile[0];
				_ue5_min = _ue_sterile[1];
				_ue6_min = _ue_sterile[2];
				_um4_min = _um_sterile[0];
				_um5_min = _um_sterile[1];
				_um6_min = _um_sterile[2];
			}
    }
  }

	std::cout << "MINS:\n"
		<< "chi2: " << _chi2_min << "\n"
		<< "m4: " << _m4_min << " m5: " << _m5_min << " m6: " << _m6_min << "\n"
		<< "ue4: " << _ue4_min << " ue5: " << _ue5_min << " ue6: " << _ue6_min << "\n"
		<< "um4: " << _um4_min << " um5: " << _um5_min << " um6: " << _um6_min << std::endl;
		

  // Create output File
  std::string outfile = rdr.tag + "_proc.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  // Create our trees
  TTree * t_app_90 = new TTree("nueapp_90","nueapp_90");
  TTree * t_app_99 = new TTree("nueapp_99","nueapp_99");
  TTree * t_m45_90 = new TTree("dm241vsdm251_90","dm241vsdm251_90");
  TTree * t_m45_99 = new TTree("dm241vsdm251_99","dm241vsdm251_99");
  TTree * t_m46_90 = new TTree("dm241vsdm261_90","dm241vsdm261_90");
  TTree * t_m46_99 = new TTree("dm241vsdm261_99","dm241vsdm261_99");


  float chi2, sin22th, dm241, dm251, dm261, chi2_min, sin22th_min, dm241_min, dm251_min, dm261_min;
  chi2_min = 99999;

  // Nue appearance
  t_app_90->Branch("chi2",&chi2,"chi2/F");
  t_app_90->Branch("sin22th",&sin22th,"sin22th/F");
  t_app_90->Branch("dm2",&dm241,"dm2/F");
  t_app_99->Branch("chi2",&chi2,"chi2/F");
  t_app_99->Branch("sin22th",&sin22th,"sin22th/F");
  t_app_99->Branch("dm2",&dm241,"dm2/F");

  float mstep = TMath::Log10(100.f/.01f)/rdr.gridpts_dm2;
  float sin22step = TMath::Log10(1.f/float(1e-5))/rdr.gridpts_sin22th;
  std::vector < std::vector < float > > chi2grid;
  chi2grid.assign(rdr.gridpts_sin22th, std::vector < float > (rdr.gridpts_dm2, 0.));

  for(int i = 0; i < v_chi2.size(); i++){
    float mysin22th = 4*pow(v_um4[i],2)*pow(v_ue4[i],2);
    int is = (int)TMath::Nint(TMath::Log10(mysin22th/1e-5)/sin22step);
    int im = (int)TMath::Nint(TMath::Log10(pow(v_mnu4[i],2)/.01)/mstep);

    if(is < rdr.gridpts_sin22th && im < rdr.gridpts_dm2 && is > -1){
      if(chi2grid[is][im] == 0)
        chi2grid[is][im] = v_chi2[i];
      else
        chi2grid[is][im] = TMath::Min(chi2grid[is][im],v_chi2[i]);
    }
    else
      continue;

    if(chi2grid[is][im] < chi2_min){
      chi2_min = chi2grid[is][im];
      dm241_min = pow(v_mnu4[i],2);
      sin22th_min = mysin22th;
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm2 min: " << dm241_min  << " sin22th min: " << sin22th_min << std::endl;

  for(int i = 0; i < rdr.gridpts_sin22th; i++){
    for(int j = 0; j < rdr.gridpts_dm2; j++){
      sin22th = pow(10,(i/float(rdr.gridpts_sin22th) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
      dm241 = pow(10,(j/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
      if(chi2grid[i][j] != 0){
        chi2 = chi2grid[i][j];

        if(chi2grid[i][j] <= chi2_min + 6.25)
          t_app_90->Fill();
        if(chi2grid[i][j] <= chi2_min + 11.34)
          t_app_99->Fill();
      }
    }
  }

  // Dm241 vs dm2 51 plot
  chi2_min = 99999999;
  t_m45_90->Branch("chi2",&chi2,"chi2/F");
  t_m45_90->Branch("dm241",&dm241,"dm241/F");
  t_m45_90->Branch("dm251",&dm251,"dm251/F");
  t_m45_99->Branch("chi2",&chi2,"chi2/F");
  t_m45_99->Branch("dm241",&dm241,"dm241/F");
  t_m45_99->Branch("dm251",&dm251,"dm251/F");

  chi2grid.assign(rdr.gridpts_dm2, std::vector < float > (rdr.gridpts_dm2, 0.));
  for(int i = 0; i < v_chi2.size(); i++){
    int im = (int)TMath::Nint(TMath::Log10(pow(v_mnu4[i],2)/.01)/mstep);
    int jm = (int)TMath::Nint(TMath::Log10(pow(v_mnu5[i],2)/.01)/mstep);

    if(im < rdr.gridpts_dm2 && jm < rdr.gridpts_dm2){
      if(chi2grid[im][jm] == 0)
        chi2grid[im][jm] = v_chi2[i];
      else
        chi2grid[im][jm] = TMath::Min(chi2grid[im][jm],v_chi2[i]);
    }
    else
      continue;
    if(chi2grid[im][jm] < chi2_min){
      chi2_min = chi2grid[im][jm];
      dm241_min = pow(v_mnu4[i],2);
      dm251_min = pow(v_mnu5[i],2);
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm241 min: " << dm241_min  << " dm251 min: " << dm251_min << std::endl;

  for(int i = 0; i < rdr.gridpts_dm2; i++){
    for(int j = 0; j < rdr.gridpts_dm2; j++){
      dm241 = pow(10,(i/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
      dm251 = pow(10,(j/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
      if(chi2grid[i][j] != 0){
        chi2 = chi2grid[i][j];

        if(chi2grid[i][j] <= chi2_min + 6.25)
          t_m45_90->Fill();
        if(chi2grid[i][j] <= chi2_min + 11.34)
          t_m45_99->Fill();
      }
    }
  }

  // Dm241 vs dm2 61 plot
  chi2_min = 99999999;
  t_m46_90->Branch("chi2",&chi2,"chi2/F");
  t_m46_90->Branch("dm241",&dm241,"dm241/F");
  t_m46_90->Branch("dm261",&dm261,"dm261/F");
  t_m46_99->Branch("chi2",&chi2,"chi2/F");
  t_m46_99->Branch("dm241",&dm241,"dm241/F");
  t_m46_99->Branch("dm261",&dm261,"dm261/F");

  chi2grid.assign(rdr.gridpts_dm2, std::vector < float > (rdr.gridpts_dm2, 0.));
  for(int i = 0; i < v_chi2.size(); i++){
    int im = (int)TMath::Nint(TMath::Log10(pow(v_mnu4[i],2)/.01)/mstep);
    int jm = (int)TMath::Nint(TMath::Log10(pow(v_mnu6[i],2)/.01)/mstep);

    if(im < rdr.gridpts_dm2 && jm < rdr.gridpts_dm2){
      if(chi2grid[im][jm] == 0)
        chi2grid[im][jm] = v_chi2[i];
      else
        chi2grid[im][jm] = TMath::Min(chi2grid[im][jm],v_chi2[i]);
    }
    else
      continue;

    if(chi2grid[im][jm] < chi2_min){
      chi2_min = chi2grid[im][jm];
      dm241_min = pow(v_mnu4[i],2);
      dm261_min = pow(v_mnu6[i],2);
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm241 min: " << dm241_min  << " dm261 min: " << dm261_min << std::endl;

  for(int i = 0; i < rdr.gridpts_dm2; i++){
    for(int j = 0; j < rdr.gridpts_dm2; j++){
      dm241 = pow(10,(i/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
      dm261 = pow(10,(j/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
			if(chi2grid[i][j] != 0){
        chi2 = chi2grid[i][j];

        if(chi2grid[i][j] <= chi2_min + 6.25)
          t_m46_90->Fill();
        if(chi2grid[i][j] <= chi2_min + 11.34)
          t_m46_99->Fill();
      }
    }
  }


  t_app_99->Write();
  t_app_90->Write();
  t_m45_90->Write();
  t_m45_99->Write();
  t_m46_90->Write();
  t_m46_99->Write();

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
