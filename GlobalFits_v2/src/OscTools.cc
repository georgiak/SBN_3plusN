#include "OscTools.h"

OutTree::OutTree(std::string tag){
  TString t_tag = tag;
  myTree = new TTree(t_tag,t_tag);
  myTree->Branch("chi2",&chi2,"chi2/F");
	myTree->Branch("dof",&dof);
	myTree->Branch("m_sterile",&m_sterile,"m_sterile[3]/F");
	myTree->Branch("um_sterile",&um_sterile,"um_sterile[3]/F");
	myTree->Branch("ue_sterile",&ue_sterile,"ue_sterile[3]/F");
	myTree->Branch("phi_sterile",&phi_sterile,"phi_sterile[3]/F");
}

void OutTree::Fill(float _chi2, float _dof, neutrinoModel _nuModel){
  chi2 = _chi2; dof = _dof;
  m_sterile[0] = _nuModel.mNu[0];    m_sterile[1] = _nuModel.mNu[1];    m_sterile[2] = _nuModel.mNu[2];
  ue_sterile[0] = _nuModel.Ue[0];    ue_sterile[1] = _nuModel.Ue[1];    ue_sterile[2] = _nuModel.Ue[2];
  um_sterile[0] = _nuModel.Um[0];    um_sterile[1] = _nuModel.Um[1];    um_sterile[2] = _nuModel.Um[2];
  phi_sterile[0] = _nuModel.phi[0];  phi_sterile[1] = _nuModel.phi[1]; 	phi_sterile[2] = _nuModel.phi[2];

  //std::cout << chi2 << " " << dof  << " " <<  m_sterile[0] << " " << ue_sterile[0] << " " <<  um_sterile[0] << std::endl;

  myTree->Fill();
}

Oscillator::Oscillator(float _dm2Min, float _dm2Max, float _UMin, float _UMax, float _USqMax, float _stepSize, float _temperature, int _nSteriles, int _gridpts, bool _CPConserving, int _nmcgen, int seed){
  dm2Min = _dm2Min;   dm2Max = _dm2Max;
  UMin = _UMin;       UMax = _UMax;
  gridpts = _gridpts;
  USqMax = _USqMax;
  nSteriles = _nSteriles;
  nMCGen = _nmcgen;
  CPConserving = _CPConserving;

  usingUe = true;
  usingUm = true;

  RanGen.SetSeed(seed);
  step = RanGen.Rndm() * _stepSize;
	temp = RanGen.Rndm() * _temperature;

  for(int i = 0; i < dm2VecMaxDim; i++){
      dm2Vec[i] = pow(10,TMath::Log10(dm2Min) + double(i) / (dm2VecMaxDim-1) * TMath::Log10(dm2Max/dm2Min));
  }
}

// Since the osc probability is something that'll be played with quite a bit (ie, when we add matter effects) let's put it over here!
oscContribution getOscContributionsNueApp(neutrinoModel model, bool nubar, bool cpv){
	oscContribution oscCon;

	oscCon.dm2[0] = pow(model.mNu[0],2);
	oscCon.aMuE[0] = 4.*model.Ue[0]*model.Um[0]*(model.Ue[0]*model.Um[0] + model.Ue[1]*model.Um[1]*cos(model.phi[0]) + model.Ue[2]*model.Um[2]*cos(model.phi[1]));
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[0] = 2.*model.Ue[0]*model.Um[0]*(model.Ue[1]*model.Um[1]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[1]));
		else    oscCon.aMuE_CPV[0] = -2.*model.Ue[0]*model.Um[0]*(model.Ue[1]*model.Um[1]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[1]));
	}

	oscCon.dm2[1] = pow(model.mNu[1],2);
	oscCon.aMuE[1] = 4.*model.Ue[1]*model.Um[1]*(model.Ue[0]*model.Um[0]*cos(model.phi[0]) + model.Ue[1]*model.Um[1] + model.Ue[2]*model.Um[2]*cos(model.phi[2]));
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[1] = 2.*model.Ue[1]*model.Um[1]*(-model.Ue[0]*model.Um[0]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[2]));
		else    oscCon.aMuE_CPV[1] = -2.*model.Ue[1]*model.Um[1]*(-model.Ue[0]*model.Um[0]*sin(model.phi[0]) + model.Ue[2]*model.Um[2]*sin(model.phi[2]));
	}

	oscCon.dm2[2] = pow(model.mNu[2],2);
	oscCon.aMuE[2] = 4.*model.Ue[2]*model.Um[2]*(model.Ue[0]*model.Um[0]*cos(model.phi[1]) + model.Ue[1]*model.Um[1]*cos(model.phi[2]) + model.Ue[2]*model.Um[2]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[2] = 2.*model.Ue[2]*model.Um[2]*(-model.Ue[0]*model.Um[0]*sin(model.phi[1]) - model.Ue[1]*model.Um[1]*sin(model.phi[2]));
		else    oscCon.aMuE_CPV[2] = -2.*model.Ue[2]*model.Um[2]*(-model.Ue[0]*model.Um[0]*sin(model.phi[1]) - model.Ue[1]*model.Um[1]*sin(model.phi[2]));
	}

	//amueboonecpv=+2.*ue(6)*um(6)*(-ue(4)*um(4)*sin(phi(2))-ue(5)*um(5)*sin(phi(3)))

	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));
	oscCon.aMuE[3] = -4.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*cos(model.phi[0]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[3] = 2.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*sin(model.phi[0]);
		else oscCon.aMuE_CPV[3] = -2.*model.Ue[0]*model.Ue[1]*model.Um[0]*model.Um[1]*sin(model.phi[0]);
	}

	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));
	oscCon.aMuE[4] = -4.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*cos(model.phi[1]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[4] = 2.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*sin(model.phi[1]);
		else    oscCon.aMuE_CPV[4] = -2.*model.Ue[0]*model.Ue[2]*model.Um[0]*model.Um[2]*sin(model.phi[1]);
	}

	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));
	oscCon.aMuE[5] = -4.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*cos(model.phi[2]);
	if(cpv) {
		if(!nubar)  oscCon.aMuE_CPV[5] = 2.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*sin(model.phi[2]);
		else    oscCon.aMuE_CPV[5] = -2.*model.Ue[1]*model.Ue[2]*model.Um[1]*model.Um[2]*sin(model.phi[2]);
	}

	return oscCon;
}

oscContribution getOscContributionsNueDis(neutrinoModel model){
	oscContribution oscCon;

	oscCon.dm2[0] = pow(model.mNu[0],2);
	oscCon.aEE[0] = -4 * pow(model.Ue[0],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[1] = pow(model.mNu[1],2);
	oscCon.aEE[1] = -4 * pow(model.Ue[1],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[2] = pow(model.mNu[2],2);
	oscCon.aEE[2] = -4 * pow(model.Ue[2],2)*(1. - pow(model.Ue[0],2) - pow(model.Ue[1],2) - pow(model.Ue[2],2));

	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));
	oscCon.aEE[3] = -4. * pow(model.Ue[0],2) * pow(model.Ue[1],2);

	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));
	oscCon.aEE[4] = -4. * pow(model.Ue[0],2) * pow(model.Ue[2],2);

	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));
	oscCon.aEE[5] = -4. * pow(model.Ue[1],2) * pow(model.Ue[2],2);

	return oscCon;
}

oscContribution getOscContributionsNumuDis(neutrinoModel model){
	oscContribution oscCon;

	oscCon.aMuMu[0] = -4. * pow(model.Um[0],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[0] = pow(model.mNu[0],2);

	oscCon.aMuMu[1] = -4. * pow(model.Um[1],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[1] = pow(model.mNu[1],2);

	oscCon.aMuMu[2] = -4. * pow(model.Um[2],2) * (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2));
	oscCon.dm2[2] = pow(model.mNu[2],2);

	oscCon.aMuMu[3] = -4. * pow(model.Um[0],2) * pow(model.Um[1],2);
	oscCon.dm2[3] = abs(pow(model.mNu[1],2) - pow(model.mNu[0],2));

	oscCon.aMuMu[4] = -4. * pow(model.Um[0],2) * pow(model.Um[2],2);
	oscCon.dm2[4] = abs(pow(model.mNu[2],2) - pow(model.mNu[0],2));

	oscCon.aMuMu[5] = -4. * pow(model.Um[2],2) * pow(model.Um[1],2);
	oscCon.dm2[5] = abs(pow(model.mNu[2],2) - pow(model.mNu[1],2));

	return oscCon;
}

// Get your models sorted out (all of this has been tested to shit. we're good.)
neutrinoModel Oscillator::InitializeMarkovParams(){

  // Initialize new model!
  neutrinoModel modelOld;
  modelOld.zero();
  bool reject;
  int nCPFactors = nSteriles*(nSteriles-1)/2;

  do{
    // Fill up our random array
    RanGen.RndmArray(13,ran);

    // Initial Params for Mass and Mixing
    if(nSteriles > 0){
      for(int i = 0; i < nSteriles; i++){
        modelOld.mNu[i] = pow(10., (TMath::Log10(dm2Min) + ran[3*i]*TMath::Log10(dm2Max/dm2Min))/2);
        //modelOld.mNu[i] = pow(10., (TMath::Log10(.317) + ran[3*i]*TMath::Log10(3.832/.316))/2);
        modelOld.Ue[i] = pow(10., TMath::Log10(UMin) + ran[3*i + 1] * TMath::Log10(UMax/UMin));
        //modelOld.Ue[i] = pow(10., TMath::Log10(.051) + ran[3*i + 1] * TMath::Log10(.57/0.051));
        modelOld.Um[i] = pow(10., TMath::Log10(UMin) + ran[3*i + 2] * TMath::Log10(UMax/UMin));
      }
    }
    // Now, let's do the CP factors, phi
    if(nCPFactors > 0){
      for(int i = 0; i < nCPFactors; i++){
        modelOld.phi[i] = double(ran[nSteriles*3+i]*2*TMath::Pi());
        if(CPConserving){
          if(modelOld.phi[i] < TMath::Pi()) modelOld.phi[i] = 0;
          else  modelOld.phi[i] = TMath::Pi();
        }
      }
    }
    reject = RejectModel(modelOld);
  }while(reject);

  return modelOld;
}

neutrinoModel Oscillator::NewModel(neutrinoModel modelOld){

  // Initialize new model!
  neutrinoModel model;
  model.zero();
  bool reject;
  int nCPFactors = nSteriles*(nSteriles-1)/2;

  do{
    // Generate some random numbers!
    RanGen.RndmArray(13,ran);

    // Alright, let's step forward with these masses and mixing matrix elements!
    for(int i = 0; i < nSteriles; i++){
			model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[3*i] - .5)*2*step*TMath::Log10(dm2Max/dm2Min))/2);
      //model.mNu[i] = pow(10., (TMath::Log10(modelOld.mNu[i]) + (ran[3*i] - .5)*2*step*TMath::Log10(3.832/.316))/2);
      //if(usingUe) model.Ue[i] = pow(10.,TMath::Log10(modelOld.Ue[i]) + (ran[3*i+1] - .5)*2*step*TMath::Log10(.57/0.051));
      if(usingUe) model.Ue[i] = pow(10.,TMath::Log10(modelOld.Ue[i]) + (ran[3*i+1] - .5)*2*step*TMath::Log10(UMax/UMin));
      else    model.Ue[i] = 0.;

      if(usingUm) model.Um[i] = pow(10.,TMath::Log10(modelOld.Um[i]) + (ran[3*i+2] - .5)*2*step*TMath::Log10(UMax/UMin));
      else    model.Um[i] = 0.;
    }
    if(nCPFactors > 0){
      for(int j = 0; j < nCPFactors; j++){
        model.phi[j] = modelOld.phi[j] + 2.*(ran[nSteriles*3 + j] - 0.5)*2.*TMath::Pi()*step;
        if(CPConserving == 1){
          if(model.phi[j] < TMath::Pi())    model.phi[j] = 0;
          else model.phi[j] = TMath::Pi();
        }
      }
    }
    reject = RejectModel(model);
		// for prospect, make sure ue4 is within proper bounds
		//reject = reject || model.Ue[0] < .051 || model.Ue[0] > .57;
  }while(reject);

  return model;
}

bool Oscillator::RejectModel(neutrinoModel model){

  int nCPFactors = nSteriles*(nSteriles-1)/2;
  // Now, we'll reject the model if matrix elements are too large
  reject1 = pow(model.Ue[0],2) + pow(model.Um[0],2) > USqMax || pow(model.Ue[1],2) + pow(model.Um[1],2) > USqMax || pow(model.Ue[2],2) + pow(model.Um[2],2) > USqMax || pow(model.Ue[0],2) + pow(model.Ue[1],2) + pow(model.Ue[2],2) > USqMax || pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2) > USqMax;

  // Another condition can be applied to avoid a negative under the square root for atmospheric neutrinos
  if(UsingAtm){
    double dmuMax = .25;

    double A = (1. - pow(model.Um[0],2) - pow(model.Um[1],2) - pow(model.Um[2],2)) * (pow(model.Um[0],2) + pow(model.Um[1],2) + pow(model.Um[2],2)) +
          pow((model.Um[0]*model.Um[1]),2) + pow((model.Um[0]*model.Um[2]),2) + pow((model.Um[1]*model.Um[2]),2);

    reject1 = reject1 || (1. - 4*A) < pow(1-2*dmuMax,2);
  }

  // More rejection conditions!
  reject2 = false;
  reject3 = false;
  reject4 = false;
  if(nSteriles > 1){
    // DEGENERACY FIX
    reject2 = model.mNu[1] < model.mNu[0] || abs(pow(model.mNu[1],2) - pow(model.mNu[0],2)) < dm2Min;
  }
  if(nSteriles > 2){
    // DEGENERACY FIX
    reject3 = model.mNu[2] < model.mNu[0] || model.mNu[2] < model.mNu[1] || abs(pow(model.mNu[2],2) - pow(model.mNu[0],2)) < dm2Min || abs(pow(model.mNu[2],2) - pow(model.mNu[1],2)) < dm2Min;
  }

  // For the Markov chain case, gotta check a few more things. Essentially whether or not we've stepped out of bounds.
  if(nSteriles > 0){
    for(int i = 0; i < nSteriles; i++){
      reject4 = reject4 || pow(model.mNu[i],2) < dm2Min || pow(model.mNu[i],2) > dm2Max;

      if(usingUe)  reject4 = reject4 || model.Ue[i] < UMin || model.Ue[i] > UMax;
      if(usingUm)  reject4 = reject4 || model.Um[i] < UMin || model.Um[i] > UMax;
    }
  }
  if(nCPFactors > 0){
    for(int i = 0; i < nCPFactors; i++){
      reject4 = reject4 || model.phi[i] < 0 || model.phi[i] > 2*TMath::Pi();
    }
  }

  return (reject1 || reject2 || reject3 || reject4);
}

double sinFunc(double x){
    // Sine function for integration
    return sin(x)/x;
}
double sineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&sinFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return ig.Integral(0.,x);
}
double cosFunc(double x){
    // cosine function for integration
    return (cos(x) - 1)/x;
}
double cosineInt(double x){
    // Sine function for integration reasons
    ROOT::Math::Functor1D wf(&cosFunc);
    ROOT::Math::Integrator ig;
    ig.SetFunction(wf);

    return TMath::EulerGamma() + log(x) + ig.Integral(0.,x);
}
