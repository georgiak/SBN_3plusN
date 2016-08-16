// This is just a temp thing to make a plot for a talk
void makeMBPlot(){

	double binso[12] = {200.,  300.,  375.,  475.,  550.,  675.,  800.,  950.,  1100.,  1300.,  1500.,  3000.};
	double bins = 11;

	TH1D* hist_obs = new TH1D("hist_obs", "", bins,binso);
	TH1D* hist_pred_31 = new TH1D("hist_pred_31", "", bins,binso);
	TH1D* hist_pred_32 = new TH1D("hist_pred_32", "", bins,binso);
	TH1D* hist_pred_33 = new TH1D("hist_pred_33", "", bins,binso);

	TCanvas* c9 = new TCanvas("c3","Signal vs Prediction",700,700);
	c9->SetLeftMargin(.15);
	c9->SetBottomMargin(.15);
	c9->SetTopMargin(.05);
	c9->SetRightMargin(.05);
	c9->cd();

	hist_obs->SetTitle("Observed vs. Predicted Signal for MiniBooNE;Energy (MeV);Events/Bin");
	hist_obs->GetXaxis()->SetTitleOffset(1);
	hist_obs->GetYaxis()->SetTitleOffset(1.4);
	hist_obs->GetXaxis()->SetTitleFont(62);
	hist_obs->GetYaxis()->SetTitleFont(62);
	hist_obs->GetYaxis()->CenterTitle();
	hist_obs->GetXaxis()->CenterTitle();
	hist_obs->GetXaxis()->SetTitleSize(0.05);
	//hist_obs->GetXaxis()->SetLabelSize(0.04);
	hist_obs->GetXaxis()->SetLabelOffset(0.001);
	hist_obs->GetYaxis()->SetTitleSize(0.05);
	//hist_obs->GetYaxis()->SetLabelSize(0.04);
	hist_obs->SetStats(kFALSE);

	globInit();
	boonePackage pack = mbNuPack;

	int nBins = 11;
	int nBins_mu = 8;
	int nFOscEvts = 17204;

	neutrinoModel threep1, threep2, threep3;
	threep1.zero();
	threep1.Ue[0] = .15; 	threep1.Um[0] = .17;	threep1.mNu[0] = sqrt(.92);
	threep2.zero();
	threep2.Ue[0] = .15; 	threep2.Um[0] = .13;	threep2.mNu[0] = sqrt(.92);
	threep2.Ue[1] = .069;	threep2.Um[1] = .16; 	threep2.mNu[1] = sqrt(17.); 	threep2.phi[0] = 1.8*TMath::Pi();
	threep3.zero();
	threep3.Ue[0] = .11; 	threep3.Um[0] = .12;	threep3.mNu[0] = sqrt(.9);
	threep3.Ue[1] = .11;	threep3.Um[1] = .17; 	threep3.mNu[1] = sqrt(17.); 	threep3.phi[0] = 1.6*TMath::Pi();
	threep3.Ue[2] = .11;	threep3.Um[2] = .14;	threep3.mNu[2] = sqrt(22.);	threep3.phi[1] = .28*TMath::Pi(); threep3.phi[2] = 1.4*TMath::Pi();

	// Initialize contributions from the oscillation probability
    oscContribution oscCont;

	_signal.resize(nBins);
	_fullData.resize(nBins + nBins_mu);
	_prediction.resize(nBins + nBins_mu);

	double minEBins[nBins];
	double maxEBins[nBins];
	double ETru[nFOscEvts], LTru[nFOscEvts];

	// THREEP1
	// Initialize vars
		for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] = 0.;
	        _fullData[iB] = pack.NueData[iB];
	        _prediction[iB] = pack.NueBgr[iB];
	    }
	    for(int iB = 0; iB < nBins_mu; iB++){
	        _fullData[iB + nBins] = pack.NumuData[iB];
	        _prediction[iB + nBins] = pack.Numu[iB];
	    }
	    full_covMatrix.ResizeTo(nBins + nBins + nBins_mu, nBins + nBins + nBins_mu);
	    full_covMatrix.Zero();
	    covMatrix.ResizeTo(nBins + nBins_mu, nBins + nBins_mu);
	    covMatrix.Zero();

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep1, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_31->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
	    }

	// threep2// Initialize vars
		for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] = 0.;
	        _fullData[iB] = pack.NueData[iB];
	        _prediction[iB] = pack.NueBgr[iB];
	    }
	    for(int iB = 0; iB < nBins_mu; iB++){
	        _fullData[iB + nBins] = pack.NumuData[iB];
	        _prediction[iB + nBins] = pack.Numu[iB];
	    }
	    full_covMatrix.ResizeTo(nBins + nBins + nBins_mu, nBins + nBins + nBins_mu);
	    full_covMatrix.Zero();
	    covMatrix.ResizeTo(nBins + nBins_mu, nBins + nBins_mu);
	    covMatrix.Zero();

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep2, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_32->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
	    }

	// threep3// Initialize vars
		for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] = 0.;
	        _fullData[iB] = pack.NueData[iB];
	        _prediction[iB] = pack.NueBgr[iB];
	    }
	    for(int iB = 0; iB < nBins_mu; iB++){
	        _fullData[iB + nBins] = pack.NumuData[iB];
	        _prediction[iB + nBins] = pack.Numu[iB];
	    }
	    full_covMatrix.ResizeTo(nBins + nBins + nBins_mu, nBins + nBins + nBins_mu);
	    full_covMatrix.Zero();
	    covMatrix.ResizeTo(nBins + nBins_mu, nBins + nBins_mu);
	    covMatrix.Zero();

	    for(int iFOsc = 0; iFOsc < nFOscEvts; iFOsc++){   // Loop over full oscillation events
	        for(int iB = 0; iB < nBins; iB++){    // Loop over energy bins to fill the prediction vector pred[]
	            minEBins[iB] = pack.EnuQE[iB];
	            maxEBins[iB] = pack.EnuQE[iB+1];

	            if(pack.FOsc_EnuQE[iFOsc] > minEBins[iB] && pack.FOsc_EnuQE[iFOsc] < maxEBins[iB]){
	                // Get prediction signal by multiplying osc prob by weight of each event

	                ETru[iFOsc] = pack.FOsc_EnuTrue[iFOsc];
	                LTru[iFOsc] = pack.FOsc_LnuTrue[iFOsc];

	                // Get oscillation probability contributions
	                oscCont = getOscContributionsNueApp(threep3, false, true);

	                for(int iContribution = 0; iContribution < 6; iContribution++){

	                    _signal[iB] += pack.FOsc_weight[iFOsc]*oscCont.aMuE[iContribution]*pow(sin(1.267*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]),2)
	                              + pack.FOsc_weight[iFOsc]*oscCont.aMuE_CPV[iContribution]*sin(1.267*2*oscCont.dm2[iContribution]*LTru[iFOsc]*.01/ETru[iFOsc]);
	                }
	            }
	        }
	    }

	    // Divide signal prediction by the number of fullosc events
	    for(int iB = 0; iB < nBins; iB++){
	        _signal[iB] /= nFOscEvts;
	    }

	    for(int iB = 0; iB < nBins; iB++){
	        _prediction[iB] += _signal[iB];
			hist_pred_33->SetBinContent(iB+1,_prediction[iB]/(binso[iB+1]-binso[iB]));
			hist_obs->SetBinContent(iB+1,_fullData[iB]/(binso[iB+1]-binso[iB]));
			hist_obs->SetBinError(iB+1,sqrt(_fullData[iB])/(binso[iB+1]-binso[iB]));
	    }

		// Okay, now draw it.
		hist_obs->SetMarkerStyle(20);
		hist_obs->SetMarkerColor(kBlack);
		hist_pred_33->SetLineWidth(3);
		hist_pred_33->SetLineColor(kMagenta);
		hist_pred_32->SetLineWidth(3);
		hist_pred_32->SetLineColor(kBlue);
		hist_pred_31->SetLineWidth(3);
		hist_pred_31->SetLineColor(kCyan);

		TLegend* legS=new TLegend(0.7,0.7,0.95,0.9);
		legS->SetFillStyle(0);
		legS->SetFillColor(0);
		legS->SetBorderSize(0);
		legS->SetTextFont(62);
		legS->SetTextSize(0.03);
		legS->AddEntry(hist_obs,"Observed","pe");
		legS->AddEntry(hist_pred_31,"3+1 Pred","f");
		legS->AddEntry(hist_pred_32,"3+2 Pred","f");
		legS->AddEntry(hist_pred_33,"3+3 Pred","f");

		hist_obs->Draw("p");
		hist_pred_31->Draw("h same");
		hist_pred_32->Draw("h same");
		hist_pred_33->Draw("h same");

		legS->Draw();

		c9->SaveAs("obsSignal.png");

}


// Old ntuplegridder implementation:

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
int mycount = 0;
int total = chi2_99->GetEntries();
for(int i0 = i0min; i0 <= i0max; i0++) for(int i1 = i1min; i1 <= i1max; i1++) for(int i2 = i2min; i2 <= i2max; i2++){

	std::cout << "Scanned " << float(mycount)/float(total) << "\% of the entries so far" << std::endl;
	std::cout << float(10000*i0 + 100*i1 + i2)/float(1000000) * 100 << "\% of the first 3dim space";
	std::cout << "This boy is " << s4Vec_99[i0][i1][i2].size() << " deep." << std::endl;
	mycount += s4Vec_99[i0][i1][i2].size();

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
