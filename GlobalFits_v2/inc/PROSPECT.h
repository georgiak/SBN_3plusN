#ifndef PROSPECT_H
#define PROSPECT_H

#include "TGraph2D.h"
#include "datasets.h"

class PROSPECT: public dataset{
  public:
    PROSPECT(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    TGraph2D * surf;
};

#endif
