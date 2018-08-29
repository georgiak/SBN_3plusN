#include "fitter.h"

using namespace std;

bool stringCompare(const string & l, const string & r){
	// We only want to compare the points, not the chisq
	return (l.substr(0,6) == r.substr(0,6));
}

int ntGridder(){

  int gridPoints = 50;
  double dmmin = 0.1;
  double dmmax = 10.;
  double umin = 0;
  double umax = .5;

	// Grab processed ntuple
  std::cout << "Loading ntuple file..." << std::endl;
	TFile *inf = new TFile("/home/dcianci/Physics/GlobalFits/SBN_3plusN/GlobalFits_v2/plotters/nt_31_all_nufact.root");
	TTree *chi2_99 = (TTree*)inf->Get("chi2_99");
	TTree *chi2_90 = (TTree*)inf->Get("chi2_90");

  float chi2[16], chi2_t, dof, step, temp, m_sterile[3], um_sterile[3], ue_sterile[3], phi_sterile[3];
  chi2_99->SetBranchAddress("chi2",&chi2);
  chi2_99->SetBranchAddress("m_sterile",&m_sterile);
  chi2_99->SetBranchAddress("um_sterile",&um_sterile);
  chi2_99->SetBranchAddress("ue_sterile",&ue_sterile);
  chi2_99->SetBranchAddress("phi_sterile",&phi_sterile);

  std::cout << chi2_99->GetEntries() << std::endl;

	// Make new file to fill with processed ntuples

	TFile *f = new TFile("nt_31_nufact_processed.root", "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

	double cl90, cl99;
  cl90 = 6.25;
  cl99 = 11.34;

  TTree *chi2_99_pr = new TTree("chi2_99","chi2_99");
  chi2_99_pr->Branch("chi2",&chi2_t,"chi2/F");
  chi2_99_pr->Branch("m_sterile",&m_sterile,"m_sterile[3]/F");
  chi2_99_pr->Branch("um_sterile",&um_sterile,"um_sterile[3]/F");
  chi2_99_pr->Branch("ue_sterile",&ue_sterile,"ue_sterile[3]/F");
  chi2_99_pr->Branch("phi_sterile",&phi_sterile,"phi_sterile[3]/F");

  TTree *chi2_90_pr = new TTree("chi2_90","chi2_90");
  chi2_90_pr->Branch("chi2",&chi2_t,"chi2[16]/F");
  chi2_90_pr->Branch("m_sterile",&m_sterile,"m_sterile[3]/F");
  chi2_90_pr->Branch("um_sterile",&um_sterile,"um_sterile[3]/F");
  chi2_90_pr->Branch("ue_sterile",&ue_sterile,"ue_sterile[3]/F");
  chi2_90_pr->Branch("phi_sterile",&phi_sterile,"phi_sterile[3]/F");

	// Now, we're going to convert to a string of coordinates.
	std::cout << "Convert entries to strings of coords" << std::endl;

	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	float ustep = (umax-umin)/float(gridPoints);
	float phistep = 2*TMath::Pi()/float(gridPoints);
	std::vector<std::string> coords;
	coords.resize(chi2_99->GetEntries());
  float chi2min = 3000;

  int _m4;	int _ue4;	int _um4;
	std::string _m4s, _m5s, _m6s, _ue4s, _ue5s, _ue6s, _um4s, _um5s, _um6s, _phi45s, _phi46s, _phi56s, _chi2s;
	for(int i = 0; i < chi2_99->GetEntries(); i++){
		chi2_99->GetEntry(i);

    // Also, we want that chi2 min
    chi2min = min(chi2[0],chi2min);

		// Now, the chi2!
		stringstream stream;
		stream << fixed << setprecision(2) << chi2[0];
    _chi2s = stream.str();

		_m4 = TMath::Max(int(floor(TMath::Log10(m_sterile[0]/dmmin)/mstep)),0);
		if(_m4 < 10) _m4s = "0" + to_string(_m4);	else _m4s = to_string(_m4);
		_ue4 = TMath::Max(int(floor(ue_sterile[0]/ustep)),0);
		if(_ue4 < 10) _ue4s = "0" + to_string(_ue4);	else _ue4s = to_string(_ue4);
		_um4 = TMath::Max(int(floor(um_sterile[0]/ustep)),0);
		if(_um4 < 10) _um4s = "0" + to_string(_um4);	else _um4s = to_string(_um4);

		coords[i] = _m4s + _ue4s + _um4s + _chi2s;
		std::cout << float(i)/float(chi2_99->GetEntries()) << "\r";
	}
	std::cout << "\nCoords all full up.\nWe've got " <<  coords.size() << " of 'em" << std::endl;
  std::cout << "Chi2min: " << chi2min << std::endl;

	sort(coords.begin(),coords.end());
	coords.erase(unique(coords.begin(), coords.end(),stringCompare), coords.end());

	std::cout << "Holy shit, that worked.\nNow, fill the new tree" << std::endl;
	for(int i = 0; i < coords.size(); i++){
		// Gotta convert those boys from string coordinates into actual positions in phase space
		_m4 = stoi(coords[i].substr(0,2));
		m_sterile[0] = pow(10,(_m4/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
    m_sterile[1] = 0;
    m_sterile[2] = 0;

		_ue4 = stoi(coords[i].substr(2,2));
		ue_sterile[0] = _ue4/float(gridPoints)*(umax-umin);
    ue_sterile[1] = 0;
    ue_sterile[2] = 0;

		_um4 = stoi(coords[i].substr(4,2));
		um_sterile[0] = _um4/float(gridPoints)*(umax-umin);
    um_sterile[1] = 0;
    um_sterile[2] = 0;

		chi2_t = stof(coords[i].substr(6));
    chi2_99_pr->Fill();
    if(chi2_t < chi2min + cl90)
      chi2_90_pr->Fill();
  }

	// Now, let's see how much this is really doing.
	std::cout << "For 99%%, we reduced " << chi2_99->GetEntries() << " events to " << chi2_99_pr->GetEntries() << std::endl;
	std::cout << "For 90%%, we reduced " << chi2_90->GetEntries() << " events to " << chi2_90_pr->GetEntries() << std::endl;

	// Save Ntuple to file
	chi2_99_pr->Write();
	chi2_90_pr->Write();
	f->Close();

  return 0;
}

int main(){
  ntGridder();
  return 1;
}
