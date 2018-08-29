#ifndef ATM_H
#define ATM_H

#include "datasets.h"

class Atm: public dataset{
  public:
    Atm(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int dmuVecMaxDim = 101;
    double *DmuVec, *DChi2Vec;
};

#endif
