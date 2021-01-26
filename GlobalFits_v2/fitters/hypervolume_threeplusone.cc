// This macro takes chi2's from selectoins of datasets, applies them to a nice grid and marginalizes  them according to sin22th(mm) or sin22th(em)

#include "fitter.h"
#include "TH3D.h"
#include <algorithm>

double IndexToValue(double _index, double _min, double _max, int _grdpts, std::string _scale="log"){
  if(_scale == "log"){
    return pow(10,log10(_min) + _index * log10(_max/_min)/(_grdpts+1));
  }
  else{
    std::cout << "Scale not yet supported" << std::endl;
    return -999;
  }
}

int ValueToIndex(double value, double _min, double _max, int _grdpts, std::string _scale="log"){
  if(_scale == "log"){
    return floor((_grdpts)*log10(value/_min)/log10(_max/_min));
  }
  else{
    std::cout << "Scale not yet supported" << std::endl;
    return -999;
  }
}

int ntupleProcess(){

  bool debug = false;

  int m41grd(80),t14grd(80),t24grd(80),t34grd(10);
  double m41max(9.0), m41min(1.f);
  double t14max(3.14/4.0), t14min(.01);
  double t24max(3.14/4.0), t24min(.01);
  double t34max(3.14/4.0), t34min(.01);

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

  // let's do nue disappearance first, so key will have two coordinates, the simplest case
  std::map<std::array<int,2>,double> chimap_2d;
  std::map<std::array<int,3>,double> chimap_3d;
  std::map<std::array<int,4>,double> chimap_4d;

  double _chi2, _m41, _theta14, _theta24, _theta34;
  int im41, it14, it24, it34;

  std::cout <<  "Alright. First, let's load up the nue disappearance chi2's (ie: 2d param space)" << std::endl;
  TFile* fProc2d = new TFile("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/Oct2020/NueDis_finer_Jan22.root","READ");
  TTree* tTotal = (TTree*)fProc2d->Get("Total");
  tTotal->SetBranchAddress("chi2",&_chi2);
  tTotal->SetBranchAddress("m41",&_m41);
  tTotal->SetBranchAddress("theta14",&_theta14);

  std::cout << "2d entries: " << tTotal->GetEntries() << std::endl;
  for(int i = 0; i < tTotal->GetEntries(); i++){
    tTotal->GetEntry(i);
    // get phase space position in generalized coordinates
    im41 = ValueToIndex(_m41,m41min,m41max,m41grd);
    it14 = ValueToIndex(_theta14,t14min,t14max,t14grd);

    if(im41 < 0 || it14 < 0 || im41 >= m41grd || it14 >= t14grd)
      continue; // get out of here if our point isn't on the grid

    std::array<int,2> hvid = {im41,it14};
    auto search = chimap_2d.find(hvid);
    if (search != chimap_2d.end()) {
      chimap_2d.at(hvid) = max(chimap_2d.at(hvid),_chi2);
    }
    else
      chimap_2d.insert({hvid,_chi2});
  }
  fProc2d->Close();

  std::cout <<  "Next, let's load up all other datasets except icecube chi2's (ie: 3d param space)" << std::endl;
  TFile* fProc3d = new TFile("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/Oct2020/AllAppPlusNumuDis_finer_Jan20.root","READ");
  tTotal = (TTree*)fProc3d->Get("Total");
  tTotal->SetBranchAddress("chi2",&_chi2);
  tTotal->SetBranchAddress("m41",&_m41);
  tTotal->SetBranchAddress("theta14",&_theta14);
  tTotal->SetBranchAddress("theta24",&_theta24);

  std::cout << "3d entries: " << tTotal->GetEntries() << std::endl;
  for(int i = 0; i < tTotal->GetEntries(); i++){
    tTotal->GetEntry(i);
    // get phase space position in generalized coordinates
    im41 = ValueToIndex(_m41,m41min,m41max,m41grd);
    it14 = ValueToIndex(_theta14,t14min,t14max,t14grd);
    it24 = ValueToIndex(_theta24,t24min,t24max,t24grd);

    if(im41 < 0 || it14 < 0 || it24 < 0 ||  im41 >= m41grd || it14 >= t14grd || it24 >= t24grd)
      continue; // get out of here if our point isn't on the grid

    std::array<int,3> hvid = {im41,it14,it24};
    auto search = chimap_3d.find(hvid);
    if (search != chimap_3d.end()) {
      chimap_3d.at(hvid) = max(chimap_3d.at(hvid),_chi2);
    }
    else
      chimap_3d.insert({hvid,_chi2});
  }
  fProc3d->Close();

  std::cout << "Okay. Now let's add the two together. We'll worry about icecube later." << std::endl;

  double chi2min(99999);
  std::array<int,4> chi2min_coords;

  //  Now we  worry about icecube.
  TFile* fProcIC = new TFile("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/Oct2020/IC_chi2grids.root","READ");
  std::vector<TH3D*> vecLogICHisto;
  for(int i = 0; i < 10; i++){
    vecLogICHisto.push_back((TH3D*)fProcIC->Get(("loghisto_"+std::to_string(i)).c_str()));
  }

  // add 'em  all together.
  double newchi2, ICchi2;
  for( auto const& x : chimap_3d ){
    std::array<int,2> hvid2d = {x.first[0],x.first[1]};
    std::array<int,3> hvid3d = {x.first[0],x.first[1],x.first[2]};

    //for(int ith34 = 0; ith34 < 1; ith34++){
    for(int ith34 = 0; ith34 < 10; ith34++){

      std::array<int,4> hvid4d = {x.first[0],x.first[1],x.first[2],ith34};
      if(m41vec.at(x.first[0]) < .1 || m41vec.at(x.first[0]) > 10.0 || t14vec.at(x.first[1]) < .01 || t14vec.at(x.first[1]) > 3.14/4 || t14vec.at(x.first[2]) < .01 || t14vec.at(x.first[2]) > 3.14/4){
        //std::cout << "Osc Params out of bounds for icecube interpolation. Returning big number" <<  std::endl;
        ICchi2 = 99999;
      }
      else
        ICchi2 = vecLogICHisto.at(ith34)->Interpolate(log10(m41vec.at(x.first[0])),log10(t14vec.at(x.first[1])),log10(t24vec.at(x.first[2])));

      double newchi2 = x.second + chimap_2d.at(hvid2d) + ICchi2;
      chimap_4d.insert({hvid4d,newchi2});
      if(newchi2 < chi2min){
        //std::cout << "NEWMIN" << std::endl;
        chi2min = newchi2;
        chi2min_coords = hvid4d;
      }
    }
  }

  std::cout << "Chi2min:" << chi2min << " at " << chi2min_coords[0] << " "  << chi2min_coords[1] << " " << chi2min_coords[2] << " " << chi2min_coords[3] << std::endl;
  std::cout << "Chi2min:" << chi2min << " at " << m41vec.at(chi2min_coords[0]) << " "  << t14vec.at(chi2min_coords[1]) << " " << t24vec.at(chi2min_coords[2]) << " " << t34vec.at(chi2min_coords[3]) << std::endl;

  // Cool. Cool. Cool.

  // Create output File
  std::string outfile = "hypervolume_1_20_2021_finer.root";
	std::cout << "Output File: " << outfile << std::endl;
	TFile *f = new TFile(outfile.c_str(), "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

  double cl90, cl99;
  cl90 = 6.25;
  cl99 = 11.34;

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


  for( auto const& x : chimap_4d ){
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

  // Let's test how we did.
  std::cout << "We started with a grid of " << chimap_4d.size() << " entries." << std::endl;
  std::cout << "for 99%, we have " << t_99CL->GetEntries() << " entries." << std::endl;
  std::cout << "for 90%, we have " << t_90CL->GetEntries() << " entries." << std::endl;

  t_99CL->Write();
  t_90CL->Write();
  f->Close();

  return 0;
}

int main(){
  ntupleProcess();
  return 0;
}
