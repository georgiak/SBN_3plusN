// This macro takes chi2's from selections of datasets and grids them according to your specifications

#include "fitter.h"
#include "TH3D.h"
#include <algorithm>

int ntupleProcess(std::string xml){

  bool debug = false;
  bool useIC = false;
  bool use3dplots = false;
  bool use2dplots = false;

  ProcessReader rdr;
  if(rdr.Load(xml))
      return 0;

  //
  //////////////////////////////////////////////////////
  // define the parameters of our grid
  int m41grd(90),t14grd(90),t24grd(90),t34grd(10);
  double m41max(10.0), m41min(.1f);
  double t14max(3.14/4.0), t14min(.01);
  double t24max(3.14/4.0), t24min(.01);
  double t34max(3.14/4.0), t34min(.01);
  std::string ic_hist_location = "/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/data/IC_chi2grids.root";
  /////////////////////////////////////////////////////

  // ok. before anything, let's sort our grid out.
  std::vector<double> m41vec, t14vec, t24vec, t34vec;
  m41vec.resize(m41grd);
  t14vec.resize(t14grd);
  t24vec.resize(t24grd);
  t34vec.resize(t34grd);

  for(int i = 0; i < m41grd; i++){
    m41vec[i] = IndexToValue(i+.5,m41min,m41max,m41grd);
  }
  for(int i = 0; i < t14grd; i++){
   t14vec[i] = IndexToValue(i+.5,t14min,t14max,t14grd);
  }
  for(int i = 0; i < t24grd; i++){
   t24vec[i] = IndexToValue(i+.5,t24min,t24max,t24grd);
  }
  for(int i = 0; i < t34grd; i++){
   t34vec[i] = IndexToValue(i+.5,t34min,t34max,t34grd);
  }


  int dof;
  std::vector < double > v_chi2, v_m41, v_theta14, v_theta24, v_theta34;

  // Initialize input tree variables
  double _chi2, _m41, _theta14, _theta24, _theta34;
  int _dof(0);

  std::cout << "Combining the chi2s from multiple datasets into one." << std::endl;
  std::cout << "If any trees are 3d, please have those listed FIRST in your XML" << std::endl;
  std::cout << "NTrees: " << rdr.GetNTrees() << std::endl;
  // This code was made to avert redundancy. Nue disappearance datasets (which are only sensitive to dm2 and theta14) need not be run across theta24

  // create vector containers by resizing to first (arbitrary) tree
  v_chi2.resize(rdr.GetTree(0)->GetEntries());
  std::fill(v_chi2.begin(), v_chi2.end(), 0);

  v_m41.resize(rdr.GetTree(0)->GetEntries());
  v_theta14.resize(rdr.GetTree(0)->GetEntries());
  v_theta24.resize(rdr.GetTree(0)->GetEntries());
  v_theta34.resize(rdr.GetTree(0)->GetEntries());

  for(int i = 0; i < rdr.GetNTrees(); i++){
    if(rdr.GetName(i) == "IceCube" && rdr.GetNTrees() != 1){  // don't do this if we're looking at ONLY icecube
      useIC = true;
      continue;
    }

    rdr.GetTree(i)->SetBranchAddress("chi2",&_chi2);
    rdr.GetTree(i)->SetBranchAddress("dof",&_dof);
    rdr.GetTree(i)->SetBranchAddress("m41",&_m41);
    rdr.GetTree(i)->SetBranchAddress("theta14",&_theta14);
    rdr.GetTree(i)->SetBranchAddress("theta24",&_theta24);
    rdr.GetTree(i)->GetEntry(10); // this does nothing important but it makes root happy.
    use3dplots = use3dplots || rdr.GetOscType(i);
    use2dplots = !(use2dplots || rdr.GetOscType(i));

    dof += _dof;

    for(int j = 0; j < rdr.GetTree(i)->GetEntries(); j++){
      rdr.GetTree(i)->GetEntry(j);

      if(rdr.GetOscType(i)>0){
        v_chi2[j] += _chi2;
        if(i==0){
          v_m41[j] = _m41;
          v_theta14[j] = _theta14;
          v_theta24[j] = _theta24;
        }
        else{
          if(v_m41[j] != _m41 || v_theta14[j] != _theta14 || v_theta24[j] != _theta24){
            std::cout << "MISALIGNMENT IN 3D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j] << " / " << _m41 << " " << v_theta14[j] << " / " << _theta14 << " " << v_theta24[j] << " / " << _theta24<< std::endl;
          }
        }
      }
      else if(use3dplots){ // if we've got a 2d case (nue dis)
        for(int x = 0; x < pow(v_chi2.size(),1.0/3.0); x++){
          v_chi2[j*100+x] += _chi2;

          if(v_m41[j*100] != _m41 || v_theta14[j*100] != _theta14){
            std::cout << "MISALIGNMENT IN 2D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j*100] << " / " << _m41 << " " << v_theta14[j*100] << " / " << _theta14 << std::endl;
          }
        }
      }
      else{ // if we only have 2d cases (nue dis)
          v_chi2[j] += _chi2;
          if(i==0){
            v_m41[j] = _m41;
            v_theta14[j] = _theta14;
            v_theta24[j] = _theta24;
          }
          else{
              if(v_m41[j] != _m41 || v_theta14[j] != _theta14){
            std::cout << "MISALIGNMENT IN 2D TREE " << i << " at  entry: " << j << std::endl;
            std::cout << v_m41[j] << " / " << _m41 << " " << v_theta14[j] << " / " << _theta14 << std::endl;
          }
        }
      }
    }
  }

  //  Now we  worry about icecube. only a special case because it's got 4 dimensions
  if(useIC) {
    std::cout << "Load icecube, the special case for it's fourth  dimension" << std::endl;
    TFile* fProcIC = new TFile(ic_hist_location.c_str(),"READ");
    std::vector<TH3D*> vecLogICHisto;
    for(int i = 0; i < 10; i++){
      vecLogICHisto.push_back((TH3D*)fProcIC->Get(("loghisto_"+std::to_string(i)).c_str()));
    }

    int size1 = v_chi2.size();

    v_chi2.resize(size1*10);
    v_m41.resize(size1*10);
    v_theta14.resize(size1*10);
    v_theta24.resize(size1*10);
    v_theta34.resize(size1*10);

    for(int i =  0; i < size1;i++){
      for(int ith34 = 1; ith34 < 10; ith34++){
        v_chi2.at(ith34*size1+i)= v_chi2.at(i)+vecLogICHisto.at(ith34)->Interpolate(log10(v_m41.at(i)),log10(v_theta14.at(i)),log10(v_theta24.at(i)));
        v_m41.at(ith34*size1+i) = v_m41.at(i);
        v_theta14.at(ith34*size1+i)= v_theta14.at(i);
        v_theta24.at(ith34*size1+i)= v_theta24.at(i);
        v_theta34.at(ith34*size1+i)= ith34;
      }
      v_chi2[i] += vecLogICHisto.at(0)->Interpolate(log10(v_m41.at(i)),log10(v_theta14.at(i)),log10(v_theta24.at(i)));
    }
    std::cout << "New size = " << v_chi2.size()  << std::endl;
  }

  // Now we've got a vector of all our points. Now let's put them into our grid!
  // Create output File
  std::string outfile = rdr.tag + "_hypervolume.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  double cl90, cl99;
  cl90 = 7.78;
  cl99 = 13.28;

  // Create our trees (one for disappearance and one for appearance)
  TTree * t_90CL = new TTree("chi2_90CL","chi2_90CL");
  TTree * t_99CL = new TTree("chi2_99CL","chi2_99CL");

  t_90CL->Branch("chi2",&_chi2);
  t_90CL->Branch("m41",&_m41);
  t_90CL->Branch("theta14",&_theta14);
  t_90CL->Branch("theta24",&_theta24);
  t_90CL->Branch("theta34",&_theta34);

  t_99CL->Branch("chi2",&_chi2);
  t_99CL->Branch("m41",&_m41);
  t_99CL->Branch("theta14",&_theta14);
  t_99CL->Branch("theta24",&_theta24);
  t_99CL->Branch("theta34",&_theta34);

  double chi2min(99999);
  if(useIC){  // 4d
    std::map<std::array<int,4>,double> chimapGrid_4d;
    // the following loop goes  through our dense chimap and respaces everything
    // taking the maximum chi2 when combining grid points
    std::cout << "Regridding points" << std::endl;
    for(int i = 0; i < v_chi2.size(); i++){
      int _x0 = ValueToIndex(v_m41.at(i),m41min,m41max,m41grd);
      int _x1 = ValueToIndex(v_theta14.at(i),t14min,t14max,t14grd);
      int _x2 = ValueToIndex(v_theta24.at(i),t24min,t24max,t24grd);
      int _x3 = v_theta34.at(i);  // theta34 is already an index
      chi2min = min(chi2min,v_chi2[i]);

      std::array<int,4> hvid4d = {_x0,_x1,_x2,_x3};
      auto search = chimapGrid_4d.find(hvid4d);
      if (search != chimapGrid_4d.end()) {
        chimapGrid_4d.at(hvid4d) = max(chimapGrid_4d.at(hvid4d),v_chi2.at(i));
      }
      else
        chimapGrid_4d.insert({hvid4d,v_chi2.at(i)});
    }
    // Print everything into a tree
    std::cout << "Printing Points" << std::endl;
    for( auto const& x : chimapGrid_4d ){
      _chi2 = x.second;
      _m41 = m41vec.at((int)x.first[0]);
      _theta14 = t14vec.at(x.first[1]);
      _theta24 = t24vec.at(x.first[2]);
      _theta34 = t34vec.at(x.first[3]);

      if( _chi2 < chi2min + cl99)
        t_99CL->Fill();
      if( _chi2 < chi2min + cl90)
        t_90CL->Fill();
    }
  }
  else{  // 3d
    std::map<std::array<int,3>,double> chimapGrid_3d;
    // the following loop goes  through our dense chimap and respaces everything
    // taking the maximum chi2 when combining grid points
    std::cout << "Regridding points" << std::endl;
    for(int i = 0; i < v_chi2.size(); i++){
      int _x0 = ValueToIndex(v_m41.at(i),m41min,m41max,m41grd);
      int _x1 = ValueToIndex(v_theta14.at(i),t14min,t14max,t14grd);
      int _x2 = ValueToIndex(v_theta24.at(i),t24min,t24max,t24grd);
      chi2min = min(chi2min,v_chi2.at(i));

      std::array<int,3> hvid3d = {_x0,_x1,_x2};
      auto search = chimapGrid_3d.find(hvid3d);
      if (search != chimapGrid_3d.end()) {
        chimapGrid_3d.at(hvid3d) = max(chimapGrid_3d.at(hvid3d),v_chi2.at(i));
      }
      else
        chimapGrid_3d.insert({hvid3d,v_chi2.at(i)});
    }
    // Print everything into a tree
    std::cout << "Printing Points" << std::endl;
    for( auto const& x : chimapGrid_3d ){
      _chi2 = x.second;
      _m41 = m41vec.at((int)x.first[0]);
      _theta14 = t14vec.at(x.first[1]);
      _theta24 = t24vec.at(x.first[2]);
      _theta34 = 0;

      if( _chi2 < chi2min + cl99)
        t_99CL->Fill();
      if( _chi2 < chi2min + cl90)
        t_90CL->Fill();
    }
  }
  // Let's test how we did.
  std::cout << "We started with a grid of " << v_chi2.size() << " entries." << std::endl;
  std::cout << "for 99%, we have " << t_99CL->GetEntries() << " entries." << std::endl;
  std::cout << "for 90%, we have " << t_90CL->GetEntries() << " entries." << std::endl;

  t_99CL->Write();
  t_90CL->Write();
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
