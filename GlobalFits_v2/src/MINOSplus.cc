// This whole thing ported from MINOS+ data release according to http://www-numi.fnal.gov/PublicInfo/forscientists.html

#include "MINOSplus.h"

int MINOSplus::Init(std::string dataLoc, bool debug){

  // Near MC RecoVtrue MINOS
  TH2D *NDNC_TrueNC_minos, *NDNC_NuMu_minos, *NDNC_BeamNue_minos, *NDNC_AppNue_minos, *NDNC_AppNuTau_minos;
  TH2D *NDCC_TrueNC_minos, *NDCC_NuMu_minos, *NDCC_BeamNue_minos, *NDCC_AppNue_minos, *NDCC_AppNuTau_minos;
  // Far MC RecoVtrue MINOS
  TH2D *FDNC_TrueNC_minos, *FDNC_NuMu_minos, *FDNC_BeamNue_minos, *FDNC_AppNue_minos, *FDNC_AppNuTau_minos;
  TH2D *FDCC_TrueNC_minos, *FDCC_NuMu_minos, *FDCC_BeamNue_minos, *FDCC_AppNue_minos, *FDCC_AppNuTau_minos;
  // MC (Unoscillated) MINOS
  TH1D *FDUnOscCC_MC_minos, *NDUnOscCC_MC_minos, *FDUnOscNC_MC_minos, *NDUnOscNC_MC_minos;
  // MC (Oscillated) MINOS
  TH1D *FDOscCC_MC_minos, *NDOscCC_MC_minos, *FDOscNC_MC_minos, *NDOscNC_MC_minos;
  // Data MINOS
  TH1D *FD_dataNC_minos, *FD_dataCC_minos, *ND_dataNC_minos, *ND_dataCC_minos;
  // Near MC RecoVtrue MINOS+
  TH2D *NDNC_TrueNC_minosPlus, *NDNC_NuMu_minosPlus, *NDNC_BeamNue_minosPlus, *NDNC_AppNue_minosPlus, *NDNC_AppNuTau_minosPlus;
  TH2D *NDCC_TrueNC_minosPlus, *NDCC_NuMu_minosPlus, *NDCC_BeamNue_minosPlus, *NDCC_AppNue_minosPlus, *NDCC_AppNuTau_minosPlus;
  // Far MC RecoVtrue MINOS+
  TH2D *FDNC_TrueNC_minosPlus, *FDNC_NuMu_minosPlus, *FDNC_BeamNue_minosPlus, *FDNC_AppNue_minosPlus, *FDNC_AppNuTau_minosPlus;
  TH2D *FDCC_TrueNC_minosPlus, *FDCC_NuMu_minosPlus, *FDCC_BeamNue_minosPlus, *FDCC_AppNue_minosPlus, *FDCC_AppNuTau_minosPlus;
  // MC (Unoscillated) MINOS+
  TH1D *FDUnOscCC_MC_minosPlus, *NDUnOscCC_MC_minosPlus, *FDUnOscNC_MC_minosPlus, *NDUnOscNC_MC_minosPlus;
  // MC (Oscillated) MINOS+
  TH1D *FDOscCC_MC_minosPlus, *NDOscCC_MC_minosPlus, *FDOscNC_MC_minosPlus, *NDOscNC_MC_minosPlus;
  // Data MINOS+
  TH1D *FD_dataNC_minosPlus, *FD_dataCC_minosPlus, *ND_dataNC_minosPlus, *ND_dataCC_minosPlus;
  // MC (Unoscillated) Joint MINOS/MINOS+
  TH1D *FDUnOscCC_MC, *NDUnOscCC_MC, *FDUnOscNC_MC, *NDUnOscNC_MC;
  // MC (Oscillated) Joint MINOS/MINOS+
  TH1D *FDOscCC_MC, *NDOscCC_MC, *FDOscNC_MC, *NDOscNC_MC;
  // Data Joint MINOS/MINOS+
  TH1D *FD_dataNC, *FD_dataCC, *ND_dataNC, *ND_dataCC;
  //Covariance Matrices -- relative variance
  TMatrixD* CoVarCC_relative;
  TMatrixD* CoVarNC_relative;
  //Covariance Matrices -- scaled, inverted
  TMatrixD* CoVarCC_inverted;
  TMatrixD* CoVarNC_inverted;
  //Beam focusing nuisance parameter weights
  TH1D *h_cc_HornIMiscal_nd_minos, *h_cc_HornIDist_nd_minos,  *h_cc_HornIMiscal_fd_minos, *h_cc_HornIDist_fd_minos;
  TH1D *h_nc_HornIMiscal_nd_minos, *h_nc_HornIDist_nd_minos,  *h_nc_HornIMiscal_fd_minos, *h_nc_HornIDist_fd_minos;
  TH1D *h_cc_HornIMiscal_nd_minosPlus, *h_cc_HornIDist_nd_minosPlus,  *h_cc_HornIMiscal_fd_minosPlus, *h_cc_HornIDist_fd_minosPlus;
  TH1D *h_nc_HornIMiscal_nd_minosPlus, *h_nc_HornIDist_nd_minosPlus,  *h_nc_HornIMiscal_fd_minosPlus, *h_nc_HornIDist_fd_minosPlus;

  // Load up data
  TFile *f = new TFile((dataLoc+"minosplus_dataRelease.root").c_str(),"READ");

  //Extract covariance matrices
  CoVarCC_relative = (TMatrixD*)f->Get("TotalCCCovar"); assert(CoVarCC_relative);
  CoVarNC_relative = (TMatrixD*)f->Get("TotalNCCovar"); assert(CoVarNC_relative);


  //Load inputs
  //Extract RecoToTrue MC simulations for MINOS
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC_minos",   NDNC_TrueNC_minos);   assert(NDNC_TrueNC_minos);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu_minos",     NDNC_NuMu_minos);     assert(NDNC_NuMu_minos);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue_minos",  NDNC_BeamNue_minos);  assert(NDNC_BeamNue_minos);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue_minos",   NDNC_AppNue_minos);   assert(NDNC_AppNue_minos);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau_minos", NDNC_AppNuTau_minos); assert(NDNC_AppNuTau_minos);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC_minos",   NDCC_TrueNC_minos);   assert(NDCC_TrueNC_minos);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu_minos",     NDCC_NuMu_minos);     assert(NDCC_NuMu_minos);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue_minos",  NDCC_BeamNue_minos);  assert(NDCC_BeamNue_minos);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue_minos",   NDCC_AppNue_minos);   assert(NDCC_AppNue_minos);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minos", NDCC_AppNuTau_minos); assert(NDCC_AppNuTau_minos);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC_minos",   FDNC_TrueNC_minos);   assert(FDNC_TrueNC_minos);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu_minos",     FDNC_NuMu_minos);     assert(FDNC_NuMu_minos);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue_minos",  FDNC_BeamNue_minos);  assert(FDNC_BeamNue_minos);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue_minos",   FDNC_AppNue_minos);   assert(FDNC_AppNue_minos);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau_minos", FDNC_AppNuTau_minos); assert(FDNC_AppNuTau_minos);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC_minos",   FDCC_TrueNC_minos);   assert(FDCC_TrueNC_minos);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu_minos",     FDCC_NuMu_minos);     assert(FDCC_NuMu_minos);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue_minos",  FDCC_BeamNue_minos);  assert(FDCC_BeamNue_minos);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue_minos",   FDCC_AppNue_minos);   assert(FDCC_AppNue_minos);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minos", FDCC_AppNuTau_minos); assert(FDCC_AppNuTau_minos);

  //Extract RecoToTrue MC simulations for MINOS+
  f->GetObject("hRecoToTrueNDNCSelectedTrueNC_minosPlus",   NDNC_TrueNC_minosPlus);   assert(NDNC_TrueNC_minosPlus);
  f->GetObject("hRecoToTrueNDNCSelectedNuMu_minosPlus",     NDNC_NuMu_minosPlus);     assert(NDNC_NuMu_minosPlus);
  f->GetObject("hRecoToTrueNDNCSelectedBeamNue_minosPlus",  NDNC_BeamNue_minosPlus);  assert(NDNC_BeamNue_minosPlus);
  f->GetObject("hRecoToTrueNDNCSelectedAppNue_minosPlus",   NDNC_AppNue_minosPlus);   assert(NDNC_AppNue_minosPlus);
  f->GetObject("hRecoToTrueNDNCSelectedAppNuTau_minosPlus", NDNC_AppNuTau_minosPlus); assert(NDNC_AppNuTau_minosPlus);

  f->GetObject("hRecoToTrueNDCCSelectedTrueNC_minosPlus",   NDCC_TrueNC_minosPlus);   assert(NDCC_TrueNC_minosPlus);
  f->GetObject("hRecoToTrueNDCCSelectedNuMu_minosPlus",     NDCC_NuMu_minosPlus);     assert(NDCC_NuMu_minosPlus);
  f->GetObject("hRecoToTrueNDCCSelectedBeamNue_minosPlus",  NDCC_BeamNue_minosPlus);  assert(NDCC_BeamNue_minosPlus);
  f->GetObject("hRecoToTrueNDCCSelectedAppNue_minosPlus",   NDCC_AppNue_minosPlus);   assert(NDCC_AppNue_minosPlus);
  f->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minosPlus", NDCC_AppNuTau_minosPlus); assert(NDCC_AppNuTau_minosPlus);

  f->GetObject("hRecoToTrueFDNCSelectedTrueNC_minosPlus",   FDNC_TrueNC_minosPlus);   assert(FDNC_TrueNC_minosPlus);
  f->GetObject("hRecoToTrueFDNCSelectedNuMu_minosPlus",     FDNC_NuMu_minosPlus);     assert(FDNC_NuMu_minosPlus);
  f->GetObject("hRecoToTrueFDNCSelectedBeamNue_minosPlus",  FDNC_BeamNue_minosPlus);  assert(FDNC_BeamNue_minosPlus);
  f->GetObject("hRecoToTrueFDNCSelectedAppNue_minosPlus",   FDNC_AppNue_minosPlus);   assert(FDNC_AppNue_minosPlus);
  f->GetObject("hRecoToTrueFDNCSelectedAppNuTau_minosPlus", FDNC_AppNuTau_minosPlus); assert(FDNC_AppNuTau_minosPlus);

  f->GetObject("hRecoToTrueFDCCSelectedTrueNC_minosPlus",   FDCC_TrueNC_minosPlus);   assert(FDCC_TrueNC_minosPlus);
  f->GetObject("hRecoToTrueFDCCSelectedNuMu_minosPlus",     FDCC_NuMu_minosPlus);     assert(FDCC_NuMu_minosPlus);
  f->GetObject("hRecoToTrueFDCCSelectedBeamNue_minosPlus",  FDCC_BeamNue_minosPlus);  assert(FDCC_BeamNue_minosPlus);
  f->GetObject("hRecoToTrueFDCCSelectedAppNue_minosPlus",   FDCC_AppNue_minosPlus);   assert(FDCC_AppNue_minosPlus);
  f->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minosPlus", FDCC_AppNuTau_minosPlus); assert(FDCC_AppNuTau_minosPlus);

  //Extract data histograms
  //MINOS
  f->GetObject("dataFDNC_minos", FD_dataNC_minos); assert(FD_dataNC_minos);
  f->GetObject("dataFDCC_minos", FD_dataCC_minos); assert(FD_dataCC_minos);

  f->GetObject("dataNDNC_minos", ND_dataNC_minos); assert(ND_dataNC_minos);
  f->GetObject("dataNDCC_minos", ND_dataCC_minos); assert(ND_dataCC_minos);

  //MINOS+
  f->GetObject("dataFDNC_minosPlus", FD_dataNC_minosPlus); assert(FD_dataNC_minosPlus);
  f->GetObject("dataFDCC_minosPlus", FD_dataCC_minosPlus); assert(FD_dataCC_minosPlus);

  f->GetObject("dataNDNC_minosPlus", ND_dataNC_minosPlus); assert(ND_dataNC_minosPlus);
  f->GetObject("dataNDCC_minosPlus", ND_dataCC_minosPlus); assert(ND_dataCC_minosPlus);


  //Construct MINOS/MINOS+ two detector data spectra
  TH1D* h2det_data_NC_minos = (TH1D*)GetTwoDetSpectrum(ND_dataNC_minos,FD_dataNC_minos);
  TH1D* h2det_data_CC_minos = (TH1D*)GetTwoDetSpectrum(ND_dataCC_minos,FD_dataCC_minos);
  TH1D* h2det_data_NC_minosPlus = (TH1D*)GetTwoDetSpectrum(ND_dataNC_minosPlus,FD_dataNC_minosPlus);
  TH1D* h2det_data_CC_minosPlus = (TH1D*)GetTwoDetSpectrum(ND_dataCC_minosPlus,FD_dataCC_minosPlus);

  //Combine MINOS & MINOS+ data spectra
  TH1D* dataCC = (TH1D*)h2det_data_CC_minos->Clone();
  dataCC->Add(h2det_data_CC_minosPlus);
  TH1D* dataNC = (TH1D*)h2det_data_NC_minos->Clone();
  dataNC->Add(h2det_data_NC_minosPlus);


  //Generate Unoscillated Spectrua
  //UnOscillated CC MC -- MINOS
  NDUnOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0,
							NDCC_TrueNC_minos,
							NDCC_NuMu_minos,
							NDCC_BeamNue_minos,
							NDCC_AppNue_minos,
							NDCC_AppNuTau_minos,
							1.04*kKmUnits);
  FDUnOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0,
							FDCC_TrueNC_minos,
							FDCC_NuMu_minos,
							FDCC_BeamNue_minos,
							FDCC_AppNue_minos,
							FDCC_AppNuTau_minos,
							735.0*kKmUnits);
  //UnOscillated NC MC -- MINOS
  NDUnOscNC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0,
							NDNC_TrueNC_minos,
							NDNC_NuMu_minos,
							NDNC_BeamNue_minos,
							NDNC_AppNue_minos,
							NDNC_AppNuTau_minos,
							1.04*kKmUnits);
  FDUnOscNC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0,
							FDNC_TrueNC_minos,
							FDNC_NuMu_minos,
							FDNC_BeamNue_minos,
							FDNC_AppNue_minos,
							FDNC_AppNuTau_minos,
							735.0*kKmUnits);
  //UnOscillated CC MC -- MINOS+
  NDUnOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0,
							NDCC_TrueNC_minosPlus,
							NDCC_NuMu_minosPlus,
							NDCC_BeamNue_minosPlus,
							NDCC_AppNue_minosPlus,
							NDCC_AppNuTau_minosPlus,
							1.04*kKmUnits);
  FDUnOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0,
							FDCC_TrueNC_minosPlus,
							FDCC_NuMu_minosPlus,
							FDCC_BeamNue_minosPlus,
							FDCC_AppNue_minosPlus,
							FDCC_AppNuTau_minosPlus,
							735.0*kKmUnits);
  //UnOscillated NC MC -- MINOS+
  NDUnOscNC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0,
							NDNC_TrueNC_minosPlus,
							NDNC_NuMu_minosPlus,
							NDNC_BeamNue_minosPlus,
							NDNC_AppNue_minosPlus,
							NDNC_AppNuTau_minosPlus,
							1.04*kKmUnits);
  FDUnOscNC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0,
							FDNC_TrueNC_minosPlus,
							FDNC_NuMu_minosPlus,
							FDNC_BeamNue_minosPlus,
							FDNC_AppNue_minosPlus,
							FDNC_AppNuTau_minosPlus,
							735.0*kKmUnits);
}

float MINOSplus::Chi2(Oscillator osc, neutrinoModel nu, bool debug){

/*
  GenerateOscillatedSpectra(my_pars);

  //Correct beam mismodeling with nuisance parameters
  ShiftBeamParameters(my_pars, f);

  //Construct MINOS/MINOS+ two detector prediction spectra
  TH1D* h2det_MC_CC_minos = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minos,FDOscCC_MC_minos);
  TH1D* h2det_MC_NC_minos = (TH1D*)GetTwoDetSpectrum(NDOscNC_MC_minos,FDOscNC_MC_minos);
  TH1D* h2det_MC_CC_minosPlus = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minosPlus,FDOscCC_MC_minosPlus);
  TH1D* h2det_MC_NC_minosPlus = (TH1D*)GetTwoDetSpectrum(NDOscNC_MC_minosPlus,FDOscNC_MC_minosPlus);

  //Combine MINOS & MINOS+ prediction spectra
  TH1D* predCC = (TH1D*)h2det_MC_CC_minos->Clone();
  predCC->Add(h2det_MC_CC_minosPlus);
  TH1D* predNC = (TH1D*)h2det_MC_NC_minos->Clone();
  predNC->Add(h2det_MC_NC_minosPlus);

  //Scale and invert covariance matrices for fitting
  CoVarCC_inverted = (TMatrixD*)ScaleCovarianceMatrix(predCC,CoVarCC_relative);
  CoVarNC_inverted = (TMatrixD*)ScaleCovarianceMatrix(predNC,CoVarNC_relative);

  Double_t chi2 = ComparePredWithData(predCC,
				      dataCC,
				      CoVarCC_inverted,
				      predNC,
				      dataNC,
				      CoVarNC_inverted,
				      my_pars.Dm232,
				      my_pars.sigma_hornImiscal_minos,
				      my_pars.sigma_hornIdist_minos,
				      my_pars.sigma_hornImiscal_minosPlus,
				      my_pars.sigma_hornIdist_minosPlus);



*/

  return 0.0f;
}




TH1D* GetTwoDetSpectrum(TH1D* hND, TH1D* hFD){

  int NDbins = hND->GetNbinsX();
  int FDbins = hFD->GetNbinsX();
  int Nbins = NDbins + FDbins;
  const int Nedges = Nbins + 1;

  Double_t edges[Nedges];

  edges[0]=0;
//////
  double shift = 40.0;	//shift bin edges of ND spectrum by maximum energy
			//of first spectrum to force increasing bin edges

  for(int i=1;i<=FDbins;i++){
    edges[i] = hFD->GetXaxis()->GetBinUpEdge(i);
  }
  for(int i=1;i<=NDbins;i++){
    edges[i+FDbins] = hND->GetXaxis()->GetBinUpEdge(i) + shift;
  }
//////
  TH1D* hSpec = new TH1D("","",Nbins,edges);
  hSpec->Sumw2();

  for(int i=1;i<=FDbins;i++){
    hSpec->SetBinContent(i,hFD->GetBinContent(i));
    hSpec->SetBinError(i,hFD->GetBinError(i));
  }
  for(int i=1;i<=NDbins;i++){
    hSpec->SetBinContent(i+FDbins,hND->GetBinContent(i));
    hSpec->SetBinError(i+FDbins,hND->GetBinError(i));
  }

  return hSpec;
}

TH1D* CreateTotalSpectrum(params my_pars,
			  TH2D* TrueNC,
			  TH2D* NuMu,
			  TH2D* BeamNue,
			  TH2D* AppNue,
			  TH2D* AppNuTau,
			  double baseline
			 )
{
  TH1D* vtruenc   = (TH1D*)CreateSpectrumComponent(my_pars, "TrueNC",   TrueNC,   baseline);
  TH1D* vnumu     = (TH1D*)CreateSpectrumComponent(my_pars, "NuMu",     NuMu,     baseline);
  TH1D* vbeamnue  = (TH1D*)CreateSpectrumComponent(my_pars, "BeamNue",  BeamNue,  baseline);
  TH1D* vappnue   = (TH1D*)CreateSpectrumComponent(my_pars, "AppNue",   AppNue,   baseline);
  TH1D* vappnutau = (TH1D*)CreateSpectrumComponent(my_pars, "AppNuTau", AppNuTau, baseline);

  TH1D* hTotal = new TH1D(*vtruenc);
  hTotal->Add(vnumu);
  hTotal->Add(vbeamnue);
  hTotal->Add(vappnue);
  hTotal->Add(vappnutau);

  return hTotal;
}

TH1D* CreateSpectrumComponent(neutrinoModel model, TString OscType, TH2D* oscDummy, Double_t baseline)
{

  TH1D* bintemplate = oscDummy->ProjectionY();
  bintemplate->Reset();

  const double k1267 = 1.26693276;

  // Loop over every true energy bin in the reco vs. true matrices, then loop over every reco energy in that bin
  // to calculate an oscillation weight for that reco energy based on the true energy.
  TAxis *Yaxis = oscDummy->GetYaxis();
  TAxis *Xaxis = oscDummy->GetXaxis();

  // Define Dm243 such that its actually Dm241 being altered.
  //41 = 43 + 32 + 21
  //43 = 41 - 32 - 21
  Double_t dm243 = 0.0;

  // SM neutrino params
  double Dm232 = 2.48524980574397732e-03,
  double Dm221 = 0.0000754,
  double th23 = 9.53091338545830835e-01,
  double th12 = 0.5540758073,
  double th13 = 0.149116,
  double deltaCP = 0.0,

  dm243 = model.dm41Sq - Dm232 - Dm221;

  for(Int_t x = 1; x <= Xaxis->GetNbins(); x++){
    Double_t OscWeight = 0.0;

    if(baseline > 0){

      // Default iterations (1 at bin center)
      Int_t n_LoverE = 1;
      Double_t LoverE[5];
      LoverE[0] = Xaxis->GetBinCenter(x);

      // This is averaging oscialltions in true energy bins - see Technical Note http://minos-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10203&version=2
      const Double_t W = Xaxis->GetBinWidth(x);
      const Double_t arg = k1267*dm243*W; // half-period of oscillation
      Double_t sample = W/2/sqrt(3);

      if(arg!=0) sample = TMath::ACos(TMath::Sin(arg)/arg)/arg*W/2;

      n_LoverE = 2;
      Double_t bc = LoverE[0]; // bin center
      LoverE[0] = bc - sample;
      LoverE[1] = bc + sample;

      const Double_t E = 1.0;

      for(int i = 0; i < n_LoverE; i++){

        // each Osctype has a different probability function
	      if(OscType == "TrueNC"){

	         OscWeight += FourFlavourNuMuToNuSProbability( E,
  						my_pars.Dm232,
							my_pars.th23,
							my_pars.Dm221,
							dm243,
							my_pars.th12,
							my_pars.th13,
							my_pars.th14,
							my_pars.th24,
							my_pars.th34,
							my_pars.deltaCP,
							0,
							my_pars.delta24,
							LoverE[i]*kKmUnits);
        }
        if(OscType == "NuMu"){

	         OscWeight += FourFlavourDisappearanceWeight( E,
  						my_pars.Dm232,
							my_pars.th23,
							my_pars.Dm221,
							dm243,
							my_pars.th12,
							my_pars.th13,
							my_pars.th14,
							my_pars.th24,
							my_pars.th34,
							my_pars.deltaCP,
							0,
							my_pars.delta24,
							LoverE[i]*kKmUnits);
        }
        if(OscType == "BeamNue"){

          OscWeight += FourFlavourNuESurvivalProbability( E,
  						my_pars.Dm232,
							my_pars.th23,
							my_pars.Dm221,
							dm243,
							my_pars.th12,
							my_pars.th13,
							my_pars.th14,
							my_pars.th24,
							my_pars.th34,
							my_pars.deltaCP,
							0,
							my_pars.delta24,
							LoverE[i]*kKmUnits);
        }
        if(OscType == "AppNue"){

          OscWeight += FourFlavourNuMuToNuEProbability( E,
  						my_pars.Dm232,
							my_pars.th23,
							my_pars.Dm221,
							dm243,
							my_pars.th12,
							my_pars.th13,
							my_pars.th14,
							my_pars.th24,
							my_pars.th34,
							my_pars.deltaCP,
							0,
							my_pars.delta24,
							LoverE[i]*kKmUnits);
        }
        if(OscType == "AppNuTau"){

          OscWeight += FourFlavourNuMuToNuTauProbability( E,
  						Dm232,
							th23,
							Dm221,
							dm243,
							th12,
							th13,
							my_pars.th14,
							my_pars.th24,
							my_pars.th34,
							deltaCP,
							0,
							my_pars.delta24,
							LoverE[i]*kKmUnits);
        }
      }
      // Now average this
      OscWeight /= n_LoverE;
    }
    else { // if baseline < 0

      if(OscType == "TrueNC")   OscWeight = 0.0;
      if(OscType == "NuMu")     OscWeight = 1.0;
      if(OscType == "BeamNue")  OscWeight = 1.0;
      if(OscType == "AppNue")   OscWeight = 0.0;
      if(OscType == "AppNuTau") OscWeight = 0.0;
    }

    // using the oscillation weight, fill a 1d histogram for each type of event with the oscillated reco energy
    for(Int_t y = 1; y <= Yaxis->GetNbins(); y++){

      Double_t sumWeights = 0;

      if(OscType == "TrueNC"){
        sumWeights += oscDummy->GetBinContent(x,y)*(1.0-OscWeight);
      }
      else{
	      sumWeights += oscDummy->GetBinContent(x,y)*(OscWeight);
      }
      Double_t currBinContents = bintemplate->GetBinContent( y );
      bintemplate->SetBinContent( y, sumWeights + currBinContents);
    }
  }
  return bintemplate;
}
