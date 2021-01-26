#include "fitter.h"

int ntupleProcess(std::string xml){

  bool debug = false;

  ProcessReader rdr;
  if(rdr.Load(xml))
    return 0;

  // Now, sum up the chi2's of everything and throw them in a vector
  float dof;
  std::vector < float > v_chi2, v_m41, v_theta14, v_theta24;

  // Initialize tree variables
  double _chi2, _m41, _theta14, _theta24, _theta34;
  int _dof;

  std::cout << "Combining the chi2s from multiple datasets into one." << std::endl;
  std::cout << "NTrees: " << rdr.GetNTrees() << std::endl;
  for(int i = 0; i < rdr.GetNTrees(); i++){
    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("dof",&_dof);
    rdr.GetTree(i)->SetBranchAddress("m41",&_m41);
    rdr.GetTree(i)->SetBranchAddress("theta14",&_theta14);
    rdr.GetTree(i)->SetBranchAddress("theta24",&_theta24);
    if(i == 0){
      v_chi2.resize(rdr.GetTree(i)->GetEntries());
      v_m41.resize(rdr.GetTree(i)->GetEntries());
      v_theta14.resize(rdr.GetTree(i)->GetEntries());
      v_theta24.resize(rdr.GetTree(i)->GetEntries());
    }
    rdr.GetTree(i)->GetEntry(10); // this does nothing important but it makes root happy.

    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);

      if(i == 0){
        if(j == 0)  dof = _dof;
        v_chi2[j] = _chi2;
        v_m41[j] = _m41;
        v_theta14[j] = _theta14;
        v_theta24[j] = _theta24;
      }
      else{
        if(j == 0)  dof += _dof;
        v_chi2[j] += _chi2;
        //if(_m41 != v_m41[0] || _theta14 != v_theta14[0]  || _theta24 != v_theta14[0]){
        //  std::cout << "BOOM " << j << " of " << rdr.GetTree(i)->GetEntries() << std::endl;
        //  std::cout << v_m41[j] << " " << _m41 << " --- " <<  v_theta14[j] << " " << v_theta14[0] << " --- " << v_theta24[j] << " " << v_theta14[0] << std::endl;
        //  std::cout << "Grids are NOT aligned and this code is not prepared to sort them out." << std::endl;
        //  return 0;
        //}
      }
    }
  }

  // Create output File
  std::string outfile = rdr.tag + "_EEproc.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  // Create our tree
  TTree * t_dis = new TTree("nuedis","nuedis");
  float chi2, sin22th, dm2, chi2_min, sin22th_min, dm2_min;
  chi2_min = 99999;

  // First, we marginalize everything into a grid of chi2s
  float mstep = TMath::Log10(100.f/.01f)/rdr.gridpts_dm2;
  float sin22step = TMath::Log10(1.f/float(1e-5))/rdr.gridpts_sin22th;
  std::vector < std::vector < float > > chi2grid;
  chi2grid.assign(rdr.gridpts_sin22th, std::vector < float > (rdr.gridpts_dm2, 0.));

  for(int i = 0; i < v_chi2.size(); i++){
    float mysin22th = pow(sin(2*v_theta14[i]),2);
    int is = (int)TMath::Nint(TMath::Log10(mysin22th/1e-5)/sin22step);
    int im = (int)TMath::Nint(TMath::Log10(pow(v_m41[i],2)/.01)/mstep);

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
      dm2_min = pow(v_m41[i],2);
      sin22th_min = mysin22th;
    }
  }

  std::cout << "chi2 min: " << chi2_min << " dm2 min: " << dm2_min  << " sin22th min: " << sin22th_min << std::endl;

  // Raster is the only option for now
  // Nueapp
  TTree * t_ras_95 = new TTree("nuedis_95","nuedis_95");
  t_ras_95->Branch("chi2",&chi2,"chi2/F");
  t_ras_95->Branch("sin22th",&sin22th,"sin22th/F");
  t_ras_95->Branch("dm2",&dm2,"dm2/F");

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
      if(chi2 > 2.706 + chi2minRaster){ //95%
      chi2 = chi2grid[sins][dm];
      //if(chi2 > 1.64 + chi2minRaster){ // 90%
        sin22th = pow(10,(sins/float(rdr.gridpts_sin22th) * TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
        dm2 = pow(10,(dm/float(rdr.gridpts_dm2) * TMath::Log10(100./.01) + TMath::Log10(.01)));
        t_ras_95->Fill();
        break;
      }
    }
  }
  t_ras_95->Write();

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
