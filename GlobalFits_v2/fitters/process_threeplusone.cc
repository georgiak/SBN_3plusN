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

    rdr.GetTree(i)->GetEntry(90);

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
    }
  }

  // Create output File
  std::string outfile = rdr.tag + "_proc.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  // Create our trees (one for disappearance and one for appearance)
  TTree * t_app_90 = new TTree("nueapp_90","nueapp_90");
  TTree * t_app_99 = new TTree("nueapp_99","nueapp_99");
  TTree * t_dis = new TTree("nuedis","nuedis");
  float chi2, sin22th, dm2, chi2_min, sin22th_min, dm2_min;
  chi2_min = 99999;

  // Numu disappearance
  t_dis->Branch("chi2",&chi2,"chi2/F");
  t_dis->Branch("sin22th",&sin22th,"sin22th/F");
  t_dis->Branch("dm2",&dm2,"dm2/F");

  for(int i = 0; i < v_chi2.size(); i++){
    chi2 = v_chi2[i];
    dm2 = pow(v_mnu[i],2);
    sin22th = 4 * pow(v_um4[i],2) * (1.f - pow(v_um4[i],2));
    t_dis->Fill();
  }
  t_dis->Write();


  // Nue disappearance
  t_app_90->Branch("chi2",&chi2,"chi2/F");
  t_app_90->Branch("sin22th",&sin22th,"sin22th/F");
  t_app_90->Branch("dm2",&dm2,"dm2/F");
  t_app_99->Branch("chi2",&chi2,"chi2/F");
  t_app_99->Branch("sin22th",&sin22th,"sin22th/F");
  t_app_99->Branch("dm2",&dm2,"dm2/F");

  float mstep = TMath::Log10(100.f/.01f)/rdr.gridpts_dm2;
  float sin22step = TMath::Log10(1.f/float(1e-5))/rdr.gridpts_sin22th;
  std::vector < std::vector < float > > chi2grid;
  chi2grid.assign(rdr.gridpts_sin22th, std::vector < float > (rdr.gridpts_dm2, 0.));

  for(int i = 0; i < v_chi2.size(); i++){
    float mysin22th = 4*pow(v_um4[i],2)*pow(v_ue4[i],2);
    int is = (int)TMath::Nint(TMath::Log10(mysin22th/1e-5)/sin22step);
    int im = (int)TMath::Nint(TMath::Log10(pow(v_mnu[i],2)/.01)/mstep);

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
      dm2_min = pow(v_mnu[i],2);
      sin22th_min = mysin22th;
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm2 min: " << dm2_min  << " sin22th min: " << sin22th_min << std::endl;

  for(int i = 0; i < rdr.gridpts_sin22th; i++){
    for(int j = 0; j < rdr.gridpts_dm2; j++){
      sin22th = pow(10,(i/float(rdr.gridpts_sin22th) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
      dm2 = pow(10,(j/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
      if(chi2grid[i][j] != 0){
        chi2 = chi2grid[i][j];

//        std::cout << "deltachi2 = " << chi2 - chi2_min << std::endl;
        if(chi2grid[i][j] <= chi2_min + 6.25)
          t_app_90->Fill();
        if(chi2grid[i][j] <= chi2_min + 11.34)
          t_app_99->Fill();
      }
    }
  }

  // If we're doing a raster scan:
  if(rdr.raster){
    TTree * t_app_95 = new TTree("nueapp_95","nueapp_95");
    t_app_95->Branch("chi2",&chi2,"chi2/F");
    t_app_95->Branch("sin22th",&sin22th,"sin22th/F");
    t_app_95->Branch("dm2",&dm2,"dm2/F");

    for(int dm = 0; dm < rdr.gridpts_dm2; dm++){
      float chi2minRaster = 3000.;
      int sinsStartRaster = 0;

      // Find minimum point for this particular dm2
      for(int sins = 0; sins < rdr.gridpts_sin22th; sins++){
        if(chi2grid[sins][dm] < chi2minRaster && chi2grid[sins][dm] >= chi2_min){
          chi2minRaster = chi2grid[sins][dm];
          sinsStartRaster = sins;
        }
      }
      for(int sins = sinsStartRaster; sins < rdr.gridpts_sin22th; sins++){
        chi2 = chi2grid[sins][dm];
        if(chi2 > 3.84 + chi2minRaster){
          sin22th = pow(10,(sins/float(rdr.gridpts_sin22th) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
          dm2 = pow(10,(dm/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
          t_app_95->Fill();
          break;
        }
      }
    }
    t_app_95->Write();
  }

  t_app_99->Write();
  t_app_90->Write();

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
