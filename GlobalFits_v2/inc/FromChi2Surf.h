#ifndef FROMCHI2SURF_H
#define FROMCHI2SURF_H

#include "datasets.h"

class FromChi2Surf: public dataset{
  public:
    FromChi2Surf(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    TH2D *mySurf;
};

#endif
