#ifndef OSCTOOLS_H
#define OSCTOOLS_H
#include "globalFit.h"

struct neutrinoModel{
    double mNu[3], Ue[3], Um[3], phi[3];
    double dm41Sq, dm51Sq, dm61Sq, dm54Sq, dm64Sq, dm65Sq;
    void zero(){
        for(int i = 0; i < 3; i ++){
            mNu[i] = 0; Ue[i] = 0;
            Um[i] = 0;  phi[i] = 0;
        }
    }
    void difference(){
        dm41Sq = pow(mNu[0],2);
        dm51Sq = pow(mNu[1],2);
        dm61Sq = pow(mNu[2],2);
        dm54Sq = dm51Sq - dm41Sq;
        dm64Sq = dm61Sq - dm41Sq;
        dm65Sq = dm61Sq - dm51Sq;
    }
};

struct oscContribution{
    double dm2[6], aMuE[6], aMuMu[6], aMuE_CPV[6], aEE[6];
};


oscContribution getOscContributionsNueApp(neutrinoModel model, bool nubar, bool cpv);
oscContribution getOscContributionsNueDis(neutrinoModel model);
oscContribution getOscContributionsNumuDis(neutrinoModel model);


class Oscillator{
  public:
    Oscillator(){};
    Oscillator(float _dm2Min, float _dm2Max, float _UMin, float _UMax, float _USqMax, float _stepSize, float _temperature, int _nSteriles, int _gridpts, bool _CPConserving, int seed);
    neutrinoModel InitializeMarkovParams();
    neutrinoModel NewModel(neutrinoModel modelOld);
    bool RejectModel(neutrinoModel model);

    int GridSize(){ return gridpts; };

  private:
    TRandom RanGen;
    float ran[13];
    float dm2Min, dm2Max, UMin, UMax, USqMax, temp, step;
    int nSteriles, gridpts;
    bool CPConserving, reject1, reject2, reject3, reject4, usingUe, usingUm;
};

class OutTree{
  public:
    OutTree(){};
    OutTree(std::string tag);

    void Fill(float _chi2, float _dof, neutrinoModel _nuModel);
    void Write(){ myTree->Write();  };
    TTree *Tree(){ return myTree->CloneTree(); };

  private:
    TTree *myTree;
    float chi2, dof, m_sterile[3], um_sterile[3], ue_sterile[3], phi_sterile[3];
};


#endif
