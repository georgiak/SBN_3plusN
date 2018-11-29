   #include "fitter.h"

int ntupleProcess(std::string xml){

  bool debug = false;
  int chainsize = 20000;

  ProcessReader rdr;
  if(rdr.Load(xml))
    return 0;


  // Just feed in the total tree
  // Initialize tree variables
  float _chi2, _m_sterile[3], _um_sterile[3], _ue_sterile[3], _dof;

  rdr.GetTree(0)->SetBranchAddress("chi2",&_chi2);
  rdr.GetTree(0)->SetBranchAddress("m_sterile",&_m_sterile);
  rdr.GetTree(0)->SetBranchAddress("um_sterile",&_um_sterile);
  rdr.GetTree(0)->SetBranchAddress("ue_sterile",&_ue_sterile);
  rdr.GetTree(0)->SetBranchAddress("dof",&_dof);

  // First, loop through everything to get all those chi2 mins
  double chi2Min = 9999999;
  double mnuMin, ue4Min, um4Min;
  std::array < double, 100 > chi2Min_chain_a, mnuMin_chain_a, ue4Min_chain_a, um4Min_chain_a;
  int ichain = -1;

  std::cout << "Loop once to find minima" << std::endl;
  for(int i = 0; i < rdr.GetTree(0)->GetEntries(); i++){
    rdr.GetTree(0)->GetEntry(i);

    // Chain local stuff
    if(i%chainsize == 0){
      ichain++;
      chi2Min_chain_a[ichain] = _chi2;
      mnuMin_chain_a[ichain] = _m_sterile[0];
      ue4Min_chain_a[ichain] = _ue_sterile[0];
      um4Min_chain_a[ichain] = _um_sterile[0];
    }
    else if(_chi2 < chi2Min_chain_a[ichain]){
      chi2Min_chain_a[ichain] = _chi2;
      mnuMin_chain_a[ichain] = _m_sterile[0];
      ue4Min_chain_a[ichain] = _ue_sterile[0];
      um4Min_chain_a[ichain] = _um_sterile[0];
    }

    if(_chi2 < chi2Min){
      chi2Min = _chi2;
      mnuMin = _m_sterile[0];
      ue4Min = _ue_sterile[0];
      um4Min = _um_sterile[0];
    }
  }

  // Now loop through again to make plots!
  // We're going to have a BUNCH of tgraphs
  std::vector < TGraph * > distLast, distMin_global, distMin_chain;
  distLast.resize(100);
  distMin_global.resize(100);
  distMin_chain.resize(100);
  ichain = -1;

  double dm_c,dm_g,dp;
  double prev_mnu = 0; double prev_ue = 0; double prev_um = 0;

  std::cout << "Loop once more to fill tgraphs" << std::endl;
  int cnt = 0;
  for(int i = 0; i < rdr.GetTree(0)->GetEntries(); i++){
    rdr.GetTree(0)->GetEntry(i);

    if(i%chainsize==0){
      ichain ++;
      distMin_chain[ichain] = new TGraph(chainsize);
      distMin_global[ichain] = new TGraph(chainsize);
      distLast[ichain] = new TGraph(chainsize);
    }

    dm_c = sqrt(pow(_m_sterile[0]-mnuMin_chain_a[ichain],2) + pow(_ue_sterile[0]-ue4Min_chain_a[ichain],2) + pow(_um_sterile[0]-um4Min_chain_a[ichain],2));
    distMin_chain[ichain]->SetPoint(cnt,i%chainsize,dm_c);

    dm_g = sqrt(pow(_m_sterile[0]-mnuMin,2) + pow(_ue_sterile[0]-ue4Min,2) + pow(_um_sterile[0]-um4Min,2));
    distMin_global[ichain]->SetPoint(cnt,i%chainsize,dm_g);

    dp = sqrt(pow(_m_sterile[0]-prev_mnu,2) + pow(_ue_sterile[0]-prev_ue,2) + pow(_um_sterile[0]-prev_um,2));
    distLast[ichain]->SetPoint(cnt,i%chainsize,dp);

    cnt++;
  }


  // Now draw three plots with allllll of these tgraphs on 'em.
  std::cout << "Draw tgraphs" << std::endl;
  TCanvas *c1 = new TCanvas("c1","Dist from chain min",1000,500);
  TCanvas *c2 = new TCanvas("c2","Dist from global min",1000,500);
  TCanvas *c3 = new TCanvas("c3","Dist from last point",1000,500);
  for(int i = 0; i < 100; i ++){
    c1->cd();
    distMin_chain[i]->Draw();
    c2->cd();
    distMin_global[i]->Draw();
    c3->cd();
    distLast[i]->Draw();
  }
  std::cout << "Print plots!" << std::endl;
  c1->SaveAs("distchain.png");
  c2->SaveAs("distglobal.png");
  c3->SaveAs("distLast.png");

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
