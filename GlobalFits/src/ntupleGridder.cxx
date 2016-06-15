/* ------------------------------------------//
Created by Davio Cianci
Jan 25th, 2016

This takes all the different markov chains and smashes them together, finds the minimum chi2 and
then outputs two final ntuples that contain all points within the 90 and 99% CL

------------------------------------------// */

#include "TLegend.h"
#include "globalFit.h"
#include "TH3D.h"
bool procOpt();

float chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56;
float m4_min,ue4_min,um4_min,m5_min,ue5_min,um5_min,m6_min,ue6_min,um6_min,phi45_min,phi46_min,phi56_min;
int steriles, nRuns, type, raster;
std::string dataset, location, output;
std::string procOptLoc;

TNtuple * chi2_99, * chi2_90, * chi2_95;
TNtuple * chi2_99_pr, * chi2_90_pr;
int gridPoints;
float dmmin, dmmax, umin, umax;
int fillProcessNT(float m4,float ue4,float um4,float m5,float ue5,float um5,float m6,float ue6,float um6,float phi45,float phi46,float phi56, int CL);

int ntGridder(){

	procOptLoc = "/Users/dcianci/Physics/SBN_3plusN/GlobalFits/inputs/";
	procOpt();

	gridPoints = 100;
	dmmin = 0.1;	dmmax = 10.;
	umin = 0;		umax = .5;

	// Grab processed ntuple
    std::cout << "Loading ntuple file..." << std::endl;
	std::string jid = Form("/nt_3%i_",steriles);
	std::string infile = output + jid + dataset + ".root";
	std::cout << "Input File: " << infile << std::endl;
	TString inputFile = infile;
	TFile *inf = new TFile(inputFile);
	TNtuple *chi2_99 = (TNtuple*)inf->Get("chi2_99");
	TNtuple *chi2_90 = (TNtuple*)inf->Get("chi2_90");


	// Make new file to fill with processed ntuples
	jid = Form("/nt_3%i_",steriles);
	std::string outfile = output + jid + dataset + "_processed.root";
	std::cout << "Output File: " << outfile << std::endl;
	TString outputFile = outfile;
	TFile *f = new TFile(outputFile, "RECREATE");
	if(f->IsZombie()){
		std::cout << "Error: couldn't create output file." << std::endl;
		return 0;
	}

	chi2_99_pr = new TNtuple("chi2_99_pr","chi2_99_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");
	chi2_90_pr = new TNtuple("chi2_90_pr","chi2_90_pr","chi2:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

	double cl90, cl99;
	if(steriles == 1){
		cl90 = 6.25;
		cl99 = 11.34;
	}
	if(steriles == 2){
		cl90 = 12.02;
		cl99 = 18.48;
	}
	if(steriles == 3){
		cl90 = 18.55;
		cl99 = 26.22;
	}

	if(steriles == 1){
		// Okay, this one's 'easy' using built-in root data structures.
		// Fill up a histogram so everything's in proper order
		std::cout << "Filling rastergram histo" << std::endl;
		TH3D * rastergram_99 = new TH3D("rg99","rg",gridPoints,.1,10.,gridPoints,0,.5,gridPoints,0,.5);
		TH3D * rastergram_90 = new TH3D("rg90","rg",gridPoints,.1,10.,gridPoints,0,.5,gridPoints,0,.5);

		// We'll deal with 99% first, yeah?
		chi2_99->SetBranchAddress("chi2",&chi2);
		chi2_99->SetBranchAddress("m4",&m4);
		chi2_99->SetBranchAddress("ue4",&ue4);
		chi2_99->SetBranchAddress("um4",&um4);

		for(int i = 0; i < chi2_99->GetEntries(); i++){
        	chi2_99->GetEntry(i);
			float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
			float ustep = (umax-umin)/float(gridPoints);
			float _m4 = ceil(TMath::Log10(m4/dmmin)/mstep);
			float _ue4 = ceil(ue4/ustep);
			float _um4 = ceil(um4/ustep);

			rastergram_99->SetBinContent(_m4,_ue4,_um4,chi2);
		}

		chi2_90->SetBranchAddress("chi2",&chi2);
		chi2_90->SetBranchAddress("m4",&m4);
		chi2_90->SetBranchAddress("ue4",&ue4);
		chi2_90->SetBranchAddress("um4",&um4);

		for(int i = 0; i < chi2_90->GetEntries(); i++){
        	chi2_90->GetEntry(i);
			float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
			float ustep = (umax-umin)/float(gridPoints);
			float _m4 = ceil(TMath::Log10(m4/dmmin)/mstep);
			float _ue4 = ceil(ue4/ustep);
			float _um4 = ceil(um4/ustep);

			rastergram_90->SetBinContent(_m4,_ue4,_um4,chi2);
		}

		// Now, fill the new ntuples
		for(int _m4 = 1; _m4 <= gridPoints; _m4++){
			for(int _ue4 = 1; _ue4 <= gridPoints; _ue4++){
				for(int _um4 = 1; _um4 <= gridPoints; _um4++){
					float m4p = pow(10,(_m4/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
					float ue4p = _ue4/float(gridPoints)*(umax-umin);
					float um4p = _um4/float(gridPoints)*(umax-umin);
					if(rastergram_99->GetBinContent(_m4,_ue4,_um4) > 0){
						chi2_99_pr->Fill(chi2,m4p,ue4p,um4p,0,0,0,0,0,0,0,0,0);
						std::cout << "m4p: " << m4p << " " << "m4: " << _m4 << std::endl;
					}
					if(rastergram_90->GetBinContent(_m4,_ue4,_um4) > 0)
						chi2_90_pr->Fill(chi2,m4p,ue4p,um4p,0,0,0,0,0,0,0,0,0);
				}
			}
		}
	}

	if(steriles == 2){
		// Okay, so this one is theoretically the hardest because we don't want to just brute-force it. So how do we go about this?
		// Cleverly is how! And just a bit at a time.

		// Okay, so initialize our 4d vectors
		std::vector <std::vector <std::vector <std::vector <int> > > > s4Vec_99;
		std::vector <std::vector <std::vector <std::vector <int> > > > s5Vec_99;
		std::vector <std::vector <int> > phiVec_99;

		// Now size them properly
		s4Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));

		// Alright, let's begin! (just 99% for now)
		chi2_99->SetBranchAddress("chi2",&chi2);
		chi2_99->SetBranchAddress("m4",&m4);
		chi2_99->SetBranchAddress("ue4",&ue4);
		chi2_99->SetBranchAddress("um4",&um4);
		chi2_99->SetBranchAddress("m5",&m5);
		chi2_99->SetBranchAddress("ue5",&ue5);
		chi2_99->SetBranchAddress("um5",&um5);
		chi2_99->SetBranchAddress("phi45",&phi45);

		// Fill up first vector (s4Vec)
		float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
		float ustep = (umax-umin)/float(gridPoints);
		float phistep = 2*TMath::Pi()/float(gridPoints);
		int i0min, i0max, i1min, i1max, i2min, i2max, i3min, i3max, i4min, i4max, i5min, i5max;		// maxes and mins for each parameter to speed shit up
		i0min = 100; i0max = 0; i1min = 100; i1max = 0; i2min = 100; i2max = 0;
		i3min = 100; i3max = 0; i4min = 100; i4max = 0; i5min = 100; i5max = 0;
		//for(int i = 0; i < chi2_99->GetEntries(); i++){
		for(int i = 550; i < 650; i++){
			chi2_99->GetEntry(i);
			int _m4 = floor(TMath::Log10(m4/dmmin)/mstep);
			int _ue4 = floor(ue4/ustep);
			int _um4 = floor(um4/ustep);

			i0min = TMath::Min(i0min,_m4);	i0max = TMath::Max(i0max,_m4);
			i1min = TMath::Min(i1min,_ue4);	i1max = TMath::Max(i1max,_ue4);
			i2min = TMath::Min(i2min,_um4);	i2max = TMath::Max(i2max,_um4);
			// Fill the appropriate spot on the grid with the entry number of our ntuple
			s4Vec_99[_m4][_ue4][_um4].push_back(i);
		}

		// Now, let's scan it.
		for(int i0 = i0min; i0 <= i0max; i0++) for(int i1 = i1min; i1 <= i1max; i1++) for(int i2 = i2min; i2 <= i2max; i2++){
			std::cout << float(10000*i0 + 100*i1 + i2)/float(1000000) * 100 << "%% \r";
			// If we have only one entry with these coordinates, fucking great! Store that away!
			if(s4Vec_99[i0][i1][i2].size() == 1){
				chi2_99->GetEntry(s4Vec_99[i0][i1][i2][0]);
				fillProcessNT(m4,ue4,um4,m5,ue5,um5,0,0,0,phi45,0,0,99);
				std::cout << "s4vec: " << m4 << " " << ue4 << " " << um4 << " " << m5 << " " << ue5 << " " << um5 << " " << phi45 << std::endl;
			}
			// Otherwise, we go deeper
			else if(s4Vec_99[i0][i1][i2].size() > 1){
				// Fill next vector
				i3min = 100; i3max = 0; i4min = 100; i4max = 0; i5min = 100; i5max = 0;
				// re-resize s5vec
				s5Vec_99.clear(); s5Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
				for(int i = 0; i < s4Vec_99[i0][i1][i2].size(); i++){
					chi2_99->GetEntry(s4Vec_99[i0][i1][i2][i]);
					int _m5 = floor(TMath::Log10(m5/dmmin)/mstep);
					int _ue5 = floor(ue5/ustep);
					int _um5 = floor(um5/ustep);

					i3min = TMath::Min(i3min,_m5);	i3max = TMath::Max(i3max,_m5);
					i4min = TMath::Min(i4min,_ue5);	i4max = TMath::Max(i4max,_ue5);
					i5min = TMath::Min(i5min,_um5);	i5max = TMath::Max(i5max,_um5);
					// Fill the appropriate spot on the grid with the entry number of our ntuple
					s5Vec_99[_m5][_ue5][_um5].push_back(s4Vec_99[i0][i1][i2][i]);
				}

				// And scan it...
				for(int i3 = i3min; i3 <= i3max; i3++) for(int i4 = i4min; i4 <= i4max; i4++) for(int i5 = i5min; i5 <= i5max; i5++){
					//If we have only one entry with these coordinates, dot dot dot
					if(s5Vec_99[i3][i4][i5].size() == 1){
						chi2_99->GetEntry(s5Vec_99[i3][i4][i5][0]);
						fillProcessNT(m4,ue4,um4,m5,ue5,um5,0,0,0,phi45,0,0,99);
						std::cout << "s5vec: " << m4 << " " << ue4 << " " << um4 << " " << m5 << " " << ue5 << " " << um5 << " " << phi45 << std::endl;

					}
					// Deeper!
					else if(s5Vec_99[i3][i4][i5].size() > 1){
						// Fill nexter vector
						phiVec_99.resize(100,vector<int>(0));
						for(int i = 0; i < s5Vec_99[i3][i4][i5].size(); i++){
							chi2_99->GetEntry(s5Vec_99[i3][i4][i5][i]);
							int _phi45 = floor(phi45/phistep);

							phiVec_99[_phi45].push_back(s5Vec_99[i3][i4][i5][i]);
						}

						for(int j = 0; j < 100; j++){
							if(phiVec_99[j].size() != 0){
								chi2_99->GetEntry(phiVec_99[j][0]);
								fillProcessNT(m4,ue4,um4,m5,ue5,um5,0,0,0,phi45,0,0,99);
								std::cout << "phi: " << m4 << " " << ue4 << " " << um4 << " " << m5 << " " << ue5 << " " << um5 << " " << phi45 << std::endl;
							}
						}
					}
				}
			}
		}
		std::cout << "99%% done!" << std::endl;

		// Now, for 90%, we don't need to do all that bullshit again. We can just go through 99% and take ones with appropriate chi2
		// find the chi2 min:
		chi2_99_pr->SetBranchAddress("chi2",&chi2);
		chi2_99_pr->SetBranchAddress("m4",&m4);
		chi2_99_pr->SetBranchAddress("ue4",&ue4);
		chi2_99_pr->SetBranchAddress("um4",&um4);
		chi2_99_pr->SetBranchAddress("m5",&m5);
		chi2_99_pr->SetBranchAddress("ue5",&ue5);
		chi2_99_pr->SetBranchAddress("um5",&um5);
		chi2_99_pr->SetBranchAddress("m6",&m6);
		chi2_99_pr->SetBranchAddress("ue6",&ue6);
		chi2_99_pr->SetBranchAddress("um6",&um6);
		chi2_99_pr->SetBranchAddress("phi45",&phi45);
		chi2_99_pr->SetBranchAddress("phi46",&phi46);
		chi2_99_pr->SetBranchAddress("phi56",&phi56);

		double chi2min = 3000;
		for(int i = 0; i < chi2_99->GetEntries(); i++){
	        chi2_99->GetEntry(i);
			if(chi2 < chi2min)
				chi2min = chi2;
	    }
		// Now fill up the 90%
		for(int i = 0; i < chi2_99_pr->GetEntries(); i++){
			chi2_99_pr->GetEntry(i);
			if(chi2 < chi2min + cl90)
				chi2_90_pr->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
		}
	}

	if(steriles == 3){
		// Okay, so this one is theoretically the hardest because we don't want to just brute-force it. So how do we go about this?
		// Cleverly is how! And just a bit at a time.

		// Okay, so initialize our 4d vectors
		std::vector <std::vector <std::vector <std::vector <int> > > > s4Vec_99;
		std::vector <std::vector <std::vector <std::vector <int> > > > s5Vec_99;
		std::vector <std::vector <std::vector <std::vector <int> > > > s6Vec_99;
		std::vector <std::vector <std::vector <std::vector <int> > > > phiVec_99;
		std::vector <std::vector <std::vector <std::vector <int> > > > s4Vec_90;
		std::vector <std::vector <std::vector <std::vector <int> > > > s5Vec_90;
		std::vector <std::vector <std::vector <std::vector <int> > > > s6Vec_90;
		std::vector <std::vector <std::vector <std::vector <int> > > > phiVec_90;
		// Now size them properly
		s4Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
		s5Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
		s6Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
		phiVec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));

		// Alright, let's begin! (just 99% for now)
		chi2_99->SetBranchAddress("chi2",&chi2);
		chi2_99->SetBranchAddress("m4",&m4);
		chi2_99->SetBranchAddress("ue4",&ue4);
		chi2_99->SetBranchAddress("um4",&um4);
		chi2_99->SetBranchAddress("m5",&m5);
		chi2_99->SetBranchAddress("ue5",&ue5);
		chi2_99->SetBranchAddress("um5",&um5);
		chi2_99->SetBranchAddress("m6",&m6);
		chi2_99->SetBranchAddress("ue6",&ue6);
		chi2_99->SetBranchAddress("um6",&um6);
		chi2_99->SetBranchAddress("phi45",&phi45);
		chi2_99->SetBranchAddress("phi46",&phi46);
		chi2_99->SetBranchAddress("phi56",&phi56);

		// Fill up first vector (s4Vec)
		float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
		float ustep = (umax-umin)/float(gridPoints);
		float phistep = 2*TMath::Pi()/float(gridPoints);
		int i0min, i0max, i1min, i1max, i2min, i2max, i3min, i3max, i4min, i4max, i5min, i5max,
		i6min, i6max, i7min, i7max, i8min, i8max, i9min, i9max, i10min, i10max, i11min, i11max;			// maxes and mins for each parameter to speed shit up
		i0min = 100; i0max = 0; i1min = 100; i1max = 0; i2min = 100; i2max = 0;
		i3min = 100; i3max = 0; i4min = 100; i4max = 0; i5min = 100; i5max = 0;
		i6min = 100; i6max = 0; i7min = 100; i7max = 0; i8min = 100; i8max = 0;
		i9min = 100; i9max = 0; i10min = 100; i10max = 0; i11min = 100; i11max = 0;
		for(int i = 0; i < chi2_99->GetEntries(); i++){
			chi2_99->GetEntry(i);
			int _m4 = floor(TMath::Log10(m4/dmmin)/mstep);
			int _ue4 = floor(ue4/ustep);
			int _um4 = floor(um4/ustep);

			i0min = TMath::Min(i0min,_m4);	i0max = TMath::Max(i0max,_m4);
			i1min = TMath::Min(i1min,_ue4);	i1max = TMath::Max(i1max,_ue4);
			i2min = TMath::Min(i2min,_um4);	i2max = TMath::Max(i2max,_um4);
			// Fill the appropriate spot on the grid with the entry number of our ntuple
			s4Vec_99[_m4][_ue4][_um4].push_back(i);
		}

		// Now, let's scan it.
		for(int i0 = i0min; i0 < i0max; i0++) for(int i1 = i1min; i1 < i1max; i1++) for(int i2 = i2min; i2 < i2max; i2++){
			std::cout << float(10000*i0 + 100*i1 + i2)/float(1000000) * 100 << "%% \r";
			// If we have only one entry with these coordinates, fucking great! Store that away!
			if(s4Vec_99[i0][i1][i2].size() == 1){
				chi2_99->GetEntry(s4Vec_99[i0][i1][i2][0]);
				fillProcessNT(m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,99);
			}
			// Otherwise, we go deeper
			if(s4Vec_99[i0][i1][i2].size() > 1){
				// Fill next vector
				i3min = 100; i3max = 0; i4min = 100; i4max = 0; i5min = 100; i5max = 0;
				s5Vec_99.clear(); s5Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
				for(int i = 0; i < s4Vec_99[i0][i1][i2].size(); i++){
					chi2_99->GetEntry(s4Vec_99[i0][i1][i2][i]);
					int _m5 = floor(TMath::Log10(m5/dmmin)/mstep);
					int _ue5 = floor(ue5/ustep);
					int _um5 = floor(um5/ustep);

					i3min = TMath::Min(i3min,_m5);	i3max = TMath::Max(i3max,_m5);
					i4min = TMath::Min(i4min,_ue5);	i4max = TMath::Max(i4max,_ue5);
					i5min = TMath::Min(i5min,_um5);	i5max = TMath::Max(i5max,_um5);
					// Fill the appropriate spot on the grid with the entry number of our ntuple
					s5Vec_99[_m5][_ue5][_um5].push_back(s4Vec_99[i0][i1][i2][i]);
				}

				// And scan it...
				for(int i3 = i3min; i3 <= i3max; i3++) for(int i4 = i4min; i4 <= i4max; i4++) for(int i5 = i5min; i5 <= i5max; i5++){
					//If we have only one entry with these coordinates, dot dot dot
					if(s5Vec_99[i3][i4][i5].size() == 1){
						chi2_99->GetEntry(s5Vec_99[i3][i4][i5][0]);
						fillProcessNT(m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,99);
					}
					// Deeper!
					if(s5Vec_99[i3][i4][i5].size() > 1){
						// Fill nexter vector
						i6min = 100; i6max = 0; i7min = 100; i7max = 0; i8min = 100; i8max = 0;
						s6Vec_99.clear(); s6Vec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
						for(int i = 0; i < s5Vec_99[i3][i4][i5].size(); i++){
							chi2_99->GetEntry(s5Vec_99[i3][i4][i5][i]);
							int _m6 = floor(TMath::Log10(m6/dmmin)/mstep);
							int _ue6 = floor(ue6/ustep);
							int _um6 = floor(um6/ustep);

							i6min = TMath::Min(i6min,_m6);	i6max = TMath::Max(i6max,_m6);
							i7min = TMath::Min(i7min,_ue6);	i7max = TMath::Max(i7max,_ue6);
							i8min = TMath::Min(i8min,_um6);	i8max = TMath::Max(i8max,_um6);
							s6Vec_99[_m6][_ue6][_um6].push_back(s5Vec_99[i3][i4][i5][0]);
						}

						// Scan it.
						for(int i6 = i6min; i6 <= i6max; i6++) for(int i7 = i7min; i7 <= i7max; i7++) for(int i8 = i8min; i8 <= i8max; i8++){
							if(s6Vec_99[i6][i7][i8].size() == 1){
								chi2_99->GetEntry(s6Vec_99[i6][i7][i8][0]);
								fillProcessNT(m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,99);
							}
							// Mas, mas
							if(s6Vec_99[i6][i7][i8].size() > 1){
								// Fill last vector
								i9min = 100; i9max = 0; i10min = 100; i10max = 0; i11min = 100; i11max = 0;
								phiVec_99.clear(); phiVec_99.resize(100,vector<vector<vector<int> > >(100,vector<vector<int>>(100)));
								for(int i = 0; i < s6Vec_99[i6][i7][i8].size(); i ++){
									chi2_99->GetEntry(s6Vec_99[i6][i7][i8][i]);
									int _phi45 = floor(phi45/phistep);
									int _phi46 = floor(phi46/phistep);
									int _phi56 = floor(phi56/phistep);

									i9min = TMath::Min(i9min,_phi45);	i9max = TMath::Max(i9max,_phi45);
									i10min = TMath::Min(i10min,_phi46);	i10max = TMath::Max(i10max,_phi46);
									i11min = TMath::Min(i11min,_phi56);	i11max = TMath::Max(i11max,_phi56);
									phiVec_99[_phi45][_phi46][_phi56].push_back(s6Vec_99[i6][i7][i8][i]);
								}

								// One last time
								for(int i9 = i9min; i9 <= i9max; i9++) for(int i10 = i10min; i10 <= i10max; i10++) for(int i11 = i11min; i11 <= i11max; i11++){
									if(phiVec_99[i9][i10][i11].size() != 0){
										chi2_99->GetEntry(phiVec_99[i9][i10][i11][0]);
										fillProcessNT(m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56,99);
									}
									// WE DID IT! Now wasn't that a doozy
								}
							}
						}
					}
				}
			}
		}
		std::cout << "99%% done!" << std::endl;

		// Now, for 90%, we don't need to do all that bullshit again. We can just go through 99% and take ones with appropriate chi2
		// find the chi2 min:
		chi2_99_pr->SetBranchAddress("chi2",&chi2);
		chi2_99_pr->SetBranchAddress("m4",&m4);
		chi2_99_pr->SetBranchAddress("ue4",&ue4);
		chi2_99_pr->SetBranchAddress("um4",&um4);
		chi2_99_pr->SetBranchAddress("m5",&m5);
		chi2_99_pr->SetBranchAddress("ue5",&ue5);
		chi2_99_pr->SetBranchAddress("um5",&um5);
		chi2_99_pr->SetBranchAddress("m6",&m6);
		chi2_99_pr->SetBranchAddress("ue6",&ue6);
		chi2_99_pr->SetBranchAddress("um6",&um6);
		chi2_99_pr->SetBranchAddress("phi45",&phi45);
		chi2_99_pr->SetBranchAddress("phi46",&phi46);
		chi2_99_pr->SetBranchAddress("phi56",&phi56);

		double chi2min = 3000;
		for(int i = 0; i < chi2_99->GetEntries(); i++){
	        chi2_99->GetEntry(i);
			if(chi2 < chi2min)
				chi2min = chi2;
	    }
		std::cout << chi2min << std::endl;
		// Now fill up the 90%
		for(int i = 0; i < chi2_99_pr->GetEntries(); i++){
			chi2_99_pr->GetEntry(i);
			if(chi2 < chi2min + cl90)
				chi2_90_pr->Fill(chi2,m4,ue4,um4,m5,ue5,um5,m6,ue6,um6,phi45,phi46,phi56);
		}
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

int fillProcessNT(float m4,float ue4,float um4,float m5,float ue5,float um5,float m6,float ue6,float um6,float phi45,float phi46,float phi56, int CL){

	// Now, we need to round to the correct place.
	float mstep = TMath::Log10(dmmax/dmmin)/float(gridPoints);
	float ustep = (umax-umin)/float(gridPoints);
	float phistep = 2*TMath::Pi()/float(gridPoints);

	int _m4 = ceil(TMath::Log10(m4/dmmin)/mstep);
	float m4p = pow(10,(_m4/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
	int _m5 = ceil(TMath::Log10(m5/dmmin)/mstep);
	float m5p = pow(10,(_m5/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));
	int _m6 = ceil(TMath::Log10(m6/dmmin)/mstep);
	float m6p = pow(10,(_m6/float(gridPoints)*TMath::Log10(dmmax/dmmin) + TMath::Log10(dmmin)));

	int _ue4 = ceil(ue4/ustep);
	float ue4p = _ue4/float(gridPoints)*(umax-umin);
	int _ue5 = ceil(ue5/ustep);
	float ue5p = _ue5/float(gridPoints)*(umax-umin);
	int _ue6 = ceil(ue6/ustep);
	float ue6p = _ue6/float(gridPoints)*(umax-umin);

	int _um4 = ceil(um4/ustep);
	float um4p = _um4/float(gridPoints)*(umax-umin);
	int _um5 = ceil(um5/ustep);
	float um5p = _um5/float(gridPoints)*(umax-umin);
	int _um6 = ceil(um6/ustep);
	float um6p = _um6/float(gridPoints)*(umax-umin);

	int _phi45 = floor(phi45/phistep);
	float phi45p = _phi45/float(gridPoints)*2*TMath::Pi();
	int _phi46 = floor(phi46/phistep);
	float phi46p = _phi46/float(gridPoints)*2*TMath::Pi();
	int _phi56 = floor(phi56/phistep);
	float phi56p = _phi56/float(gridPoints)*2*TMath::Pi();

	if(CL == 90){
		chi2_90_pr->Fill(chi2,m4p,ue4p,um4p,m5p,ue5p,um5p,m6p,ue6p,um6p,phi45p,phi46p,phi56p);
		return 1;
	}
	if(CL == 99){
		chi2_99_pr->Fill(chi2,m4p,ue4p,um4p,m5p,ue5p,um5p,m6p,ue6p,um6p,phi45p,phi46p,phi56p);
		return 1;
	}

	return 0;
}

bool procOpt(){
    // Here, we're going to read out the procOpt.txt file and assign those parameters.

    // Fill up paraVal vector
    std::string line;
    ifstream file;
    file.open(procOptLoc+"processOptions.txt");

	std::string key;
	std::string value;

	std::getline(file,line);
	std::istringstream is_line1(line);
	std::getline(is_line1,key,'=');
	std::getline(is_line1,value);
	steriles = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line2(line);
	std::getline(is_line2,key,'=');
	std::getline(is_line2,value);
	nRuns = atoi(value.c_str());

	std::getline(file,line);
	std::istringstream is_line3(line);
	std::getline(is_line3,key,'=');
	std::getline(is_line3,dataset);

	std::getline(file,line);
	std::istringstream is_line4(line);
	std::getline(is_line4,key,'=');
	std::getline(is_line4,location);

	std::getline(file,line);
	std::istringstream is_line5(line);
	std::getline(is_line5,key,'=');
	std::getline(is_line5,output);

    return true;
}

#ifndef __CINT__
int main()
{
    ntGridder();
    return 0;
}
# endif

void ntupleGridder(){
    ntGridder();
    return;
}
